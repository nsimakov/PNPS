//
// C++ Interface: pmfcalculation
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PMFCALCULATION_H
#define PMFCALCULATION_H
#include <vector>

#ifdef HARLEM_MOD
#  include <string>
#  include "mapio.h"
#endif

class TiXmlElement;
class ContWorld;
class BuildWorldNI;
class PoissonSolver;
class GOAtoms;
class VectorIntField3D;

/**
	@author Nikolay Simakov <nsimakov@andrew.cmu.edu>
*/
class PMFCalculation
{
	friend class PoissonSolver;
	public:
		PMFCalculation();
		~PMFCalculation();
		virtual int InitZero();//!<initialize zero values
		virtual int Clear();//!< delete all internal arrays and objects
	protected:
		ContWorld* m_ContWorld;
		ContWorld* m_ContWorldRough;
		//BuildWorldNI* m_BuildWorldNI;
		//PoissonSolver* m_PoissonSolver;
		
		enum{PMFPntsTypeGrid=1,PMFPntsTypePDB=3,PMFPntsTypeXYZ=2};
		//int ** PntsToCalc;
		//int * NPntsToCalc;
		int NIons;
		int *Npnts;
		GOAtoms** IonsForPMF;
		int PMFPntsType;
		int **PMFPnts;//!<Points To Calculate, in grid style
		int *PMFPntsStart;//!<Position from which start calculate (in case if calculation was stop and restart is needed)
		float **PMFxyz[3];
		double** PMF;//!<Calculated PMF
		bool** PMFCalc;//!<true if point was calculated
		
		bool   GproteinLoad;
		int    GproteinQnum;
		float *GproteinQst;
		float *GproteinPhi;
		int   *GproteinQndx;
		int LoadGprotein(const char * filename);
		bool   GproteinRoughLoad;
		int    GproteinRoughQnum;
		float *GproteinRoughQst;
		float *GproteinRoughPhi;
		int   *GproteinRoughQndx;
		int LoadGproteinRough(const char * filename);
		bool   GionLoad;
		int    *GionQnum;
		float **GionQst;
		float **GionPhi;
		int   **GionQndx;
		int LoadGion(const char * filename);
		int *grdpntQst;
		double *Gsubtract;
		
		bool WriteMapsForEachPosition;
		
		char TempResultsOut[256];
		char TempResultsOutRough[256];
	public:
		//int RunFromFile(const char *FileName);
		int Run(const TiXmlElement *RootElt);
		int RunOld(const TiXmlElement *RootElt);
		int RunSimple(const TiXmlElement *RootElt);
		int RunNoIonSizeNoProteinCharge(const TiXmlElement *RootElt);
		int RunCashing(const TiXmlElement *RootElt);
		int RunSimpleFocusing(const TiXmlElement *RootElt);
	private:
		static BuildWorldNI* CmdBuildWorld(const TiXmlElement *Elt,ContWorld* _ContWorld);
		//! Load Calculated Points
		int LoadCalculatedPoints(const char *FileName);
		int LoadCalculatedPointsOnePntPerNode(const char *FileName);
		int RearrangePoints();
		int ReDistrebutePntsOverProcesses();
		int GetPntPos(int ion,int node);
		//! Simple stright forvard
		int SimplePMFCalculation(const TiXmlElement *PoissonElt);
		int SimplePMFCalculationMPI0(const TiXmlElement *PMFElt, const TiXmlElement *PoissonElt);
		
		int BuildContWorld(ContWorld* world,BuildWorldNI* _BuildWorldNI);
		
		//!Calculations with savint intermediate results
		//int CachingPMFCalculation(const TiXmlElement *PoissonElt);
		//int CachingPMFCalculationMPI(const TiXmlElement *PMFElt,const TiXmlElement *PoissonElt);
		//int CachingPMFCalculationMPI0(const TiXmlElement *PMFElt,const TiXmlElement *PoissonElt);
		
		int PreBuildContWorld(ContWorld* world,BuildWorldNI* _BuildWorldNI);
		int FinalizePreBuildContWorld(ContWorld* world,BuildWorldNI* _BuildWorldNI,GOAtoms* goIon);
		
		//!Print result of single calculation
		void PrintSIPResultSingle(int ion, int iPos,double E,double dE,double Conv);
		void PrintSIPResultSingleRough(int ion, int iPos,double E,double dE,double Conv);
	protected:
		float** m_Surf;
		int** m_iVtmp;
		float* m_Qstat;
	public:
		//! with intermidiate arrays dropping
		//int CachingPMFCalculationMPI1(const TiXmlElement *PMFElt,const TiXmlElement *PoissonElt);
		int PreBuildContWorld1(ContWorld* world,BuildWorldNI* _BuildWorldNI);
		int FinalizePreBuildContWorld1(ContWorld* world,BuildWorldNI* _BuildWorldNI,GOAtoms* goIon);
		int IdSurf[3];
		int IdiVtmp[3];
		int IdQstat;
	protected:
	public:
		int LoadIonsForPMF(const TiXmlElement *Elt);
		int LoadPointsForPMF(const TiXmlElement *Elt);
		int ReadPointsFromPDB(int ion,const char *filename);
		int SetIonCoor(int IonType,int iPos);
		int SetIonCoor(GOAtoms* goIon, int IonType,int iPos);
		//int ProcessEnergy();
	protected:
		//!Message For MPI
		enum
		{
			MES_WHO=111,
			MES_WHATTODO=112,
			MES_WANT_CALC=113,
			MES_HAVE_RESULTS=114,
			MES_DONE=115,
			MES_HAVETOCALC=116,
			MES_ION=117,
			MES_PNT=118,
			MES_PMF=119
		};
		
};
extern PMFCalculation * PMFCalculation0;
class VectorField3D;

/**This Class is design for manipulation with PMF
 */
class PMFProcessing
{
	public:
		enum{MaskPMFCalc=1,MaskPMFInterpol=2,MaskPMFfromPoisson=3,MaskRegionBoarder=4,MaskPMFInterpolTemp=5};
		PMFProcessing();
		~PMFProcessing();
		virtual int InitZero();//!<initialize zero values
		virtual int Clear();//!< delete all internal arrays and objects
		int SetContWorld(ContWorld* _ContWorld);
		int CmdPMFProcess(const TiXmlElement *Elt);
		
		//!Read energy from file with 5 columns
#ifdef HARLEM_MOD
		
#endif
		int ReadEnergyFrom5ColumnFile(const char * FileName,double** dPMF,int** Mask,float ConvFac);
		int ReadEnergyFrom9ColumnFile(const char * FileName,double** dPMF,int** Mask,float ConvFac);
		int ReadEnergyFrom9ColumnFileCopyOctants(const char * FileName,double** dPMF,int** Mask,float ConvFac);
		int ReadEnergyFrom9ColumnFileSubEion(const char * FileName,double** dPMF,int** Mask,float ConvFac,double** Eion);
		int SetEionPolynom4(double** dPMF,int** Mask,double a4,double a3,double a2,double a1,double a0,double Rlim);
		int SetEionPolynom6(double** dPMF,int** Mask,double a6,double a5,double a4,double a3,double a2,double a1,double a0,double Rlim);
		int SetEionPolynom8(double** dPMF,int** Mask,double a8,double a7,double a6,double a5,double a4,double a3,double a2,double a1,double a0,double Rlim);
		int CalculateAvrEion(bool IonsHaveSameSize);
		int CombineEprotionEprotEion(bool IonsHaveSameSize);
	public:
		ContWorld* m_ContWorld;
		int Nions;
		int GridSize[3];
		int GS_X,GS_Y,GS_Z,GS_XY,GS_XYZ;
		double GridScale;
		
		int SetGrid(int gsX,int gsY,int gsZ,float gridscale,int m_Nion);
		
		double **PMF;//!< for PMF and Prot-Ion temporal storage
		double **Eion;//!< for Ion Energy
		double Eprot;//!<for Protein Energy
		int **Mask;//!<field for identification where the point was taken
		int **MaskIon;//!<field for identification where the point was taken
		double *EavrIon;
		
		//void SetMask()
		
		int CmdPMFProcessOld(const TiXmlElement *Elt);
		bool CmdPMFProcess_CheckPnt(int pnt);
		int CmdPMFProcess_Interpolate0(double** V,int **_Mask);
		int CmdPMFProcess_Interpolate0_1(int shft,int i,int j,double** V,int **_Mask);
		
		int CmdPMFProcess_Interpolate1(double** V,int **_Mask);
		
		int CmdConvertPoissonPotToPMF(double** V,int **_Mask,int MaskValueToSet,bool onlyWhereDnonZero);
		
		int ShowProperties();
		
		int SetNbhoodToCalcNode(double* V,int *_Mask,int Rstep);
		int SetNbhoodToCalcNodeAllIons(double** _V,int **_Mask,int Rstep);
		int LineInterpolation(double* V,int *_Mask);
		int LineInterpolation0AllIons(double** _V,int **_Mask,int Rstep);
		int LineInterpolation0(double* V,int *_Mask,int Rstep);
		int LineInterpolation1AllIons(double** _V,int **_Mask,int Rstep);
		int LineInterpolation1(double** V,int **_Mask,int ion,int Rstep);
		int LineInterpolationWithDiff(double** V,int **_Mask,int ion,int Rstep);
		
		int AverageThrLaplasAllIons(double** V,int **_Mask,int iter,int MaskNotToDo);
		int AverageThrLaplas(double* V,int *_Mask,int iter,int MaskNotToDo);
		int AverageThrLaplas1AllIons(double** V,int **_Mask,int iter,int MaskAmong);
		int AverageThrLaplas1(double* V,int *_Mask,int iter,int MaskAmong);
		int AverageThrLaplasDAllIons(double** V,int **_Mask,int iter,int MaskAmong);
		int AverageThrLaplasDnoZeroAllIons(double** V,int **_Mask,int iter,int MaskAmong);
		int SetPMFatDzero(double** V,int **_Mask,int ion,double PMFatDzero,int MaskValueToSet,int MaskWhereSet);
		int SetPMFatDnonzero(double** V,int **_Mask,int ion,double PMFatDzero,int MaskValueToSet,int MaskWhereSet);
		int SetPMFatInd(double** V,int **_Mask,int ion,double PMFVal, int MaskWhereSet);
		
		int SetSIPwhereDnotZero(double** V,int **_Mask,int ion,double SIP,float Z0,float Z1,int MaskValueToSet,int MaskWhereSet);
		
		int ConvertIndex(int **_Mask,int MaskFrom,int MaskTo);
		int ConvertIndexWithinIntZ(int **_Mask,int MaskFrom,int MaskTo,int iZ0,int iZ1);
		
		int SetBorder(int **_Mask,int Mask0,int MaskWith,int MaskValToSet);
		
		int AnalizePMF(double** V,int **_Mask);
		int RemovePMFwithHighdPMF(double** V,int **_Mask,double dPMF,int MaskToSet,int cycles);
		int RemoveLargePMFfromPNP(double** V,int **_Mask,double LargePMF,int MaskToSet);
		int RemoveBadDiffusionPoints();
		int RemoveDiffusionPointsAtHighC(float HighC);
		int CorrectFlexDiffusionWithNIndexDiff();
		
		int RemoveLargePMF(double** V,int **_Mask,int ion,double PMFVal);
		
		int CopyPMFtoContWorld();
		
		int LineInterpolation1_VF(VectorField3D *VF_PMF,VectorIntField3D *VF_Mask,int ion,int Rstep);
		int SetNbhoodToCalcNode_VF(VectorField3D *VF_PMF,VectorIntField3D *VF_Mask,int ion,int Rstep);
		int SetMaskToWherePMF_VF(VectorField3D *VF_PMF,float val,VectorIntField3D *VF_Mask,int MaskToSet);
		int SetPMFatDzero_VF(VectorField3D *VF_PMF,float val,VectorIntField3D *VF_Mask,int ion,float PMFatDzero,int MaskValueToSet,int MaskWhereSet);
};
class PMFProcTools
{
	public:
		enum{MaskPMFCalc=1,MaskPMFInterpol=2,MaskPMFfromPoisson=3,MaskRegionBoarder=4,MaskPMFInterpolTemp=5};
		PMFProcTools();
		~PMFProcTools();
		static int SetMaskWherePMFNotEqVal(VectorIntField3D* Mask,int SetMask,VectorField3D* PMF,float NotEqVal);
		static int AddToMaskWherePMFNotEqVal(VectorIntField3D* Mask,int AddMask,VectorField3D* PMF,float NotEqVal);
		static int SmoothLevelOffThrLaplas(VectorField3D* PMF,VectorIntField3D* Mask,int iters,int iMaskWhereSmooth);
};
//!Class which look for points which is interesting to calculate
class SIPPointsSearch
{
	public:
		SIPPointsSearch();
		~SIPPointsSearch();
		virtual int InitZero();//!<initialize zero values
		virtual int Clear();//!< delete all internal arrays and objects
	public:
		ContWorld* m_ContWorld;
		int SetContWorld(ContWorld* _ContWorld);
		int CmdSIPPointsSearch(const TiXmlElement *Elt);
		int ChannelPoints(int ion, float fx,float fy,float fz1,float fz2,float frstep,const char *FileName);
		int AlongZ(float fx,float fy,float fz1,float fz2,float frstep,const char *FileName);
		int DZeroBoarder(float fx,float fy,float fz1,float fz2,float fr,int ion,const char *FileName);
		
		int Box(float fx0,float fy0,float fz0, float fx1,float fy1,float fz1, float frstep,int ion,const char *FileName);
		int Cylinder(float fx,float fy,float fz1,float fz2,float fr,float frstep,int ion,const char *FileName);
				
		int CombineFiles(const char *FileName0, const char *FileName1,const char *FileNameOut,float dr);
		
};
bool r0withinr1(float x0, float y0, float z0,float x1, float y1, float z1,float dr);
#endif
