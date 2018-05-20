//
// C++ Interface: pnpsapp
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PNPSAPP_H
#define PNPSAPP_H


class TiXmlElement;
class ContWorld;
class DFTFWorld;
class PoissonSolver;
class PoissonSolverOnCuda;
class PoissonSolverOnCudaDouble;
class NernstPlankSolver;
class PoissonNernstPlanckSolver;
class BuildWorldNI;
class VectorField3D;
class PoissonBoltzmannSolver;
class PBwithLJSolver;
class IAVCalc;
class SISSGausPS;

/**
	@author Nikolay Simakov <nsimakov@andrew.cmu.edu>
*/
class SingleWorldCalcController
{
	public:
		SingleWorldCalcController();
		~SingleWorldCalcController();
		
		friend class PMFCalculation;//!<PMFCalculation use this class to build system
		friend class PoissonNernstPlanckMultiGridSolver;
		friend class IAVCalc;
		virtual int InitZero();//!<initialize zero values
		virtual int Clear();//!< delete all internal arrays and objects
	public:
		ContWorld* m_ContWorld;
		PoissonSolver* m_PoissonSolver;
		NernstPlankSolver* m_NernstPlankSolver;
		PoissonNernstPlanckSolver* m_PoissonNernstPlanckSolver;
		PoissonBoltzmannSolver* m_PoissonBoltzmannSolver;
		PBwithLJSolver* m_PBwithLJSolver;
		
		//DFTFWorld* m_DFTFWorld;
	public:
		int RunFromFile(const char *FileName);
		int RunFromCstr(const char *cstr);
		int Run(const TiXmlElement *RootElt);
	public:
		
		//int RunSingle(const TiXmlElement *RootElt);
		//int RunMPI(const TiXmlElement *RootElt);
		static ContWorld* CmdContWorld(const TiXmlElement *Elt);
		static ContWorld* CmdContWorldPoissonFocusing(const TiXmlElement *EltContWorld,const TiXmlElement *EltContWorldRough,const TiXmlElement *EltBuildWorld,const TiXmlElement *EltPoissonRough);
		static int CmdBuildWorld(const TiXmlElement *Elt,ContWorld* _ContWorld);
		static int CmdBuildWorldEu(const TiXmlElement *Elt,ContWorld* _ContWorld);
		static int CmdBuildWorldCmp(const TiXmlElement *Elt,ContWorld* _ContWorld);
		static int CmdBuildWorldScaled(const TiXmlElement *Elt,ContWorld* _ContWorld);
		static PoissonSolver* CmdPoissonSolver(const TiXmlElement *Elt,ContWorld* _ContWorld);
		static PoissonSolverOnCuda* CmdPoissonSolverOnCuda(const TiXmlElement *Elt,ContWorld* _ContWorld);
		static PoissonSolverOnCudaDouble* CmdPoissonSolverOnCudaDouble(const TiXmlElement *Elt,ContWorld* _ContWorld);
		static SISSGausPS* CmdSISSGausPS(const TiXmlElement *Elt,ContWorld* _ContWorld);
		static PoissonBoltzmannSolver* CmdPoissonBoltzmannSolver(const TiXmlElement *Elt,ContWorld* _ContWorld);
		static PBwithLJSolver* CmdPBwithLJSolver(const TiXmlElement *Elt,ContWorld* _ContWorld);
		static NernstPlankSolver* CmdNernstPlankSolver(const TiXmlElement *Elt,ContWorld* _ContWorld);
		//static int CmdPMFCalculation(const TiXmlElement *Elt);
		static int CmdMapIO(const TiXmlElement *Elt,ContWorld* _ContWorld);
		static int CmdSetDZeroWithIndex(ContWorld* _ContWorld);
		static int CmdSetDNIZeroWithDFlex(ContWorld* _ContWorld);
		static int CmdSetCatNI(ContWorld* _ContWorld,float NewC);
		static int CmdPotentialOnBoarderToZero(ContWorld* _ContWorld);
		
		int CmdPMFProcess(const TiXmlElement *Elt);
		
		static int CmdCalcRMSD(const TiXmlElement *Elt,ContWorld* _ContWorld);
		static int CmdCalcRMSDfromFile(const char *filename,ContWorld* _ContWorld,bool CompWithBC);
		static int CmdSubtractPotential(const TiXmlElement *Elt,ContWorld* _ContWorld);
		static PoissonNernstPlanckSolver* CmdPoissonNernstPlanckSolver(const TiXmlElement *Elt,ContWorld* _ContWorld);
		static int CmdAddPotential(const TiXmlElement *Elt,ContWorld* _ContWorld);
		static int CmdAddDiffPotentialToDiffGroup(const TiXmlElement *Elt,ContWorld* _ContWorld);
		static int CmdApplyAsymConc(const TiXmlElement *Elt,ContWorld* _ContWorld);
		static int CmdApplyAsymConcByScaling(const TiXmlElement *Elt,ContWorld* _ContWorld);
};

#endif
