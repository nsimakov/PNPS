//
// C++ Implementation: poissonnernstplancksolver
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#include "poissonnernstplancksolver.h"
#include "nernstplanksolver.h"
#include "poissonsolver.h"
#include "tinyxml.h"
#include "pnpdebug.h"
#include "contworld.h"
#include "math.h"
#include "pnpconstants.h"
#include "mapio.h"
#include "calcctrl.h"
#include "pnputil.h"

#include <Python.h>

PoissonNernstPlanckSolver::PoissonNernstPlanckSolver()
{
	PMFWeightModeStr.push_back("None");
	PMFWeightModeStr.push_back("Linear");
	InitZero();
}


PoissonNernstPlanckSolver::~PoissonNernstPlanckSolver()
{
	
	Clear();
}
int PoissonNernstPlanckSolver::InitZero()
{
	HaObject::SetName("PoissonNernstPlanckSolver");
	World=NULL;
	Poisson=NULL;
	NernstPlank=NULL;
	ItotErr=5550.0;
	IposErr=5550.0;
	InegErr=5550.0;
	positiveCurrentProfile=NULL;
	negativeCurrentProfile=NULL;
	currentDimension=2;
	ConvergenceCheck=10;PMFWeightMode=0;
	bLimitCurrentCalc=false;
	LimitCurrentCalcZ[0]=0;
	LimitCurrentCalcZ[1]=0;
	SaveMemory=false;
	bDouble=false;
	verbose=true;
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckSolver::Clear()
{
	DeleteCArray(positiveCurrentProfile);
	DeleteCArray(negativeCurrentProfile);
	DeleteObjByPnt(Poisson);
	DeleteObjByPnt(NernstPlank);
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckSolver::SaveXML(TiXmlElement* Elt, HaContext* p_ctxt )
{
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckSolver::LoadXML(const TiXmlElement* Elt, HaContext* p_ctxt )
{
	if(strncmp(HaObject::GetCStrName(),Elt->Value(),25))
	{
		fprintf(stderr,"ERROR: Wrong XML Element %s, expecting %s\n",Elt->Value(),HaObject::GetCStrName());
		return EXIT_FAILURE;
	}
	Clear();
	int i,gridPoint;
	//Read Primary Parameters
	Elt->GetIntAttribute("MaxIterations",&MaxIterations);
	if(Elt->GetBoolAttribute("Verbose",&verbose)!=EXIT_SUCCESS)verbose=true;
	if(Elt->GetIntAttribute("ConvergenceCheck",&ConvergenceCheck)!=EXIT_SUCCESS)ConvergenceCheck=10;
	if(Elt->GetStdStrIndex("PMFWeightMode",&PMFWeightMode,PMFWeightModeStr)!=EXIT_SUCCESS)PMFWeightMode=0;
	if(Elt->GetBoolAttribute("SaveMemory",&SaveMemory)!=EXIT_SUCCESS)SaveMemory=false;
	currentDimension=2;
	//Scale Parameters to Internal Units
	//Calculate Derivative Parameters
	if(Poisson==NULL)
	{
		Poisson = new PoissonSolver();
		Poisson->LoadXML(Elt->FirstChildElement("PoissonSolver"));
	}
	else Poisson->LoadXML(Elt->FirstChildElement("PoissonSolver"));
  
	if(NernstPlank==NULL)
	{
		NernstPlank = new NernstPlankSolver();
		NernstPlank->LoadXML(Elt->FirstChildElement("NernstPlankSolver"));
	}
	else NernstPlank->LoadXML(Elt->FirstChildElement("NernstPlankSolver"));
	
	if(Elt->GetArrOfFloatAttribute("LimitCurrentCalcZ",LimitCurrentCalcZ,2)==EXIT_SUCCESS)
	{
		bLimitCurrentCalc=true;
	}
	else
	{
		bLimitCurrentCalc=false;
		LimitCurrentCalcZ[0]=0;
		LimitCurrentCalcZ[1]=0;
	}
	
	if(Elt->GetBoolAttribute("bDouble",&bDouble)!=EXIT_SUCCESS)
		bDouble=false;
	
	Poisson->verbose=verbose;
	NernstPlank->verbose=verbose;
	
	
	ShowParameters();
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckSolver::LoadParamFromPyDict(PyObject *dict)
{
	Clear();
	int i,gridPoint;
	//Read Primary Parameters
	// =haPyDict_GetItemValueAsBool(dict,"",);
	// =haPyDict_GetItemValueAsInt(dict,"",);
	// =haPyDict_GetItemValueAsFloat(dict,"",);

	MaxIterations=haPyDict_GetItemValueAsInt(dict,"MaxIterations",MaxIterations);
	verbose=haPyDict_GetItemValueAsBool(dict,"Verbose",verbose);
	ConvergenceCheck=haPyDict_GetItemValueAsInt(dict,"ConvergenceCheck",ConvergenceCheck);
	PMFWeightMode=haPyDict_GetItemValueAsInt(dict,"PMFWeightMode",PMFWeightMode);
	SaveMemory=haPyDict_GetItemValueAsBool(dict,"SaveMemory",SaveMemory);

	currentDimension=2;
	//Scale Parameters to Internal Units
	//Calculate Derivative Parameters
	if(Poisson==NULL)
	{
		Poisson = new PoissonSolver();
		Poisson->LoadParamFromPyDict(PyDict_GetItemString(dict,"P_Param"));
	}
	else Poisson->LoadParamFromPyDict(PyDict_GetItemString(dict,"P_Param"));
  
	if(NernstPlank==NULL)
	{
		NernstPlank = new NernstPlankSolver();
		NernstPlank->LoadParamFromPyDict(PyDict_GetItemString(dict,"NP_Param"));
	}
	else NernstPlank->LoadParamFromPyDict(PyDict_GetItemString(dict,"NP_Param"));
	
	PyObject *p=NULL;
	if((p=PyDict_GetItemString(dict,"LimitCurrentCalcZ"))!=NULL)
	{
		bLimitCurrentCalc=true;
		haPy_SetCFloatArrFromListOfFloat(p,LimitCurrentCalcZ);
	}
	else
	{
		bLimitCurrentCalc=false;
		LimitCurrentCalcZ[0]=0;
		LimitCurrentCalcZ[1]=0;
	}
	bDouble=haPyDict_GetItemValueAsBool(dict,"bDouble",bDouble);
	
	Poisson->verbose=verbose;
	NernstPlank->verbose=verbose;
	
	
	ShowParameters();
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckSolver::SetContWorld(ContWorld* _world)
{
	int i;
	World=_world;
	if(Poisson!=NULL)
	{
		Poisson->SetContWorld(World);
	}
	if(NernstPlank!=NULL)
	{
		NernstPlank->SetContWorld(World);
	}
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckSolver::InitSolver()
{
	pnpPrint("\nPoisson-Nernst-Plank Solver\n");
	int iteration;
	float Tolerance;
	int i,j,k,i1,i2,gridPoint;
	float gridScale;
	int * gridSize;
  
	float *PotentialTMP, *Potential;
	float om1,om2;
	int iCurrentDimension;
	float currentSum,currentCoordinate;
	float currentDev,current;
	float currentDevPos,currentPos;
	float currentDevNeg,currentNeg;
	float maxCurrent,maxCurrentCoordinate;
	
	iCurrentDimension = currentDimension;

	gridScale = World->GridScale;
	gridSize = World->GridSize;
	Potential = World->Potential;
	
	int IDem;
	if( World->MyRank==0) IDem=World->GridSizeGlobal[iCurrentDimension]-1;
	else IDem=gridSize[iCurrentDimension]-1;
	
	if(!(positiveCurrentProfile = new double[IDem])){
		pnpError("ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(negativeCurrentProfile = new double[IDem])){
		pnpError("ERROR 204: No memory available\n");
		exit(204);
	}
	if(World->PMF!=NULL)pnpPrint("Calculate with PMF correction\n");
	
	if(World->D!=NULL&&World->NIndexing!=NULL)
		World->NIndexing->SetIonAccess(World->D);
	
	if(World->C!=NULL&&World->NIndexing!=NULL)Poisson->SetQmobForPNP();
	
	NernstPlank->InitSolver();
	Poisson->InitSolver();
	//if(SaveMemory)
	//	DeleteObjByPnt(World->NIndexing);
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckSolver::Solve()
{
	if(bDouble)
		SolveDouble();
	else
		SolveSingle();
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckSolver::SolveSingle()
{
	pnpPrint("\nPoisson-Nernst-Plank Solver\n");
	int iteration;
	float Tolerance;
	int i,j,k,i1,i2,gridPoint;
	float gridScale;
	int * gridSize;
	float *PotentialTMP, *Potential;
  //assert(pnpSolverData!=NULL);
#ifdef MPI_PARALLEL
	MPI::Status  status;
#endif
	
	gridScale = World->GridScale;
	gridSize = World->GridSize;
  //PotentialTMP = World->PotentialTMP;
	Potential = World->Potential;
	
	
	for(iteration=1;iteration<=MaxIterations;iteration++) 
	{
		int PoissonStatus,NernstPlankStatus;

		assert(Poisson);
		PoissonStatus=Poisson->Solve();
		if(PoissonStatus!=EXIT_SUCCESS)
		{
			pnpPrintGroup0("PNPSR -------------------------------------------------------------------------\n");
			pnpPrintGroup0("PNPSR  ERROR: Poisson solver has failed.\n");
			pnpPrintGroup0("PNPSR =========================================================================\n");
			return EXIT_FAILURE;
		}
		
		if(PMFWeightMode==1)
		{
			NernstPlank->PMFWeight=double(iteration)/double(MaxIterations);
			pnpPrint0("PMFWeight=%f\n",NernstPlank->PMFWeight);
		}
		NernstPlankStatus=NernstPlank->Solve();
    //NernstPlankStatus=NernstPlank->NernstPlanckSolverN2tmp();
		if(Poisson->solver!=Poisson->ArrayDirect)
			Poisson->SetQmobFromConcentration();
		if(NernstPlankStatus!=EXIT_SUCCESS)
		{
			pnpPrintGroup0("PNPSR -------------------------------------------------------------------------\n");
			pnpPrintGroup0("PNPSR  ERROR: Nernst-Planck solver has failed.\n");
			pnpPrintGroup0("PNPSR =========================================================================\n");
			return EXIT_FAILURE;
		}


		if(World->MyRank==0)
		{
			if(iteration==1)
			{
				pnpPrintGroup0("PNPSR =========================================================================\n");
				pnpPrintGroup0("PNPSR  %5s %16s %22s %12s %12s\n","Iter.", "MaxConcChange, M", "Energy,kT","dE","rel.E");
				pnpPrintGroup0("PNPSR -------------------------------------------------------------------------\n");
			}
			if(iteration!=1 && iteration%ConvergenceCheck==1)
			{
				pnpPrintGroup0("PNPSR -------------------------------------------------------------------------\n");
				pnpPrintGroup0("PNPSR  %5s %16s %22s %12s %12s\n","Iter.", "MaxConcChange, M", "Energy,kT","dE","rel.E");
				pnpPrintGroup0("PNPSR -------------------------------------------------------------------------\n");
			}
			pnpPrintGroup0("PNPSR  %5d %16.4e %22.14e %12.4e %12.4e\n", iteration, NernstPlank->MaxChange, Poisson->totalEnergy, Poisson->totalChange, Poisson->relativeChange);
			//pnpPrint("<PNPiteration Nit=\"%d\" MaxChange=\"%16.8g\" E=\"%16.8lg\" dE=\"%16.8lg\" rel_dE=\"%16.8lg\"/>\n", iteration, NernstPlank->MaxChange, Poisson->totalEnergy, Poisson->totalChange, Poisson->relativeChange);
		}
		if(iteration%ConvergenceCheck==0)
		{
			CombineAndPrintCurrents(iteration);
		}
		if(!pnpsapp->HaveTimeToRun())
		{
			pnpPrint("Run out of time\n");
			break;
		}
	}
	if(iteration%ConvergenceCheck!=0)
		CombineAndPrintCurrents(iteration);

	pnpPrintGroup0("PNPSR -------------------------------------------------------------------------\n");
	pnpPrintGroup0("PNPSR  Poisson-Nernst-Planck solver reached max iterations.\n");
	pnpPrintGroup0("PNPSR =========================================================================\n");
  //VectorField3D* J=NernstPlank->CalcJ();
  //static int NumJ=0;
  //char fout[32];
  //sprintf(fout,"jmap%d.gz",NumJ);
  //J->WriteToFile(fout);
  //NumJ++;
  //delete J;
	return EXIT_SUCCESS;
}
VectorField3D* PoissonNernstPlanckSolver::CalcCartI(int ion)
{
	pnpPrint("PoissonNernstPlanckSolver::CalcCartI ion=%d\n",ion);
	int * GS;
	int GS_X,GS_Y,GS_Z;
	int GS_XY;
	int GS_XYZ;
	int gridPoint;
	float gridScale;
	float fpoh;
	float flcon;
	float ftot;
	
	int i,j,k;
	int gridp,gridn;
	float pot;
	int gridInc;
	float *C[2];
	float * potential;
	float *UTMP=NULL;
	float *TMP=NULL;
	float **D=World->D;
	int kgrid,jgrid;
	
	NodeIndex* NIndex=World->NIndexing->NIndex;
	NodeIndexing* NIndexing=World->NIndexing;
	NodeIndex DiffBoarderMask=NodeIndexing::DiffIon0BoarderMask;
	NodeIndex BlackAndWhiteMask=NodeIndexing::BlackAndWhiteMask;
	
	GS = World->GridSize;
	GS_X = GS[0];
	GS_Y = GS[1];
	GS_Z = GS[2];
	GS_XY = GS[0]*GS[1];
	GS_XYZ = GS[2]*GS_XY;
	gridScale = World->GridScale;
	
	C[0] = World->C[0];
	C[1] = World->C[1];
	potential = World->Potential;
	
	if(NernstPlank->UTMPSingle==NULL)
		NernstPlank->UTMPSingle=new float[GS_XYZ];
	if(NernstPlank->TMPSingle==NULL)
		NernstPlank->TMPSingle=new float[GS_XYZ];
	
	UTMP=NernstPlank->UTMPSingle;
	TMP=NernstPlank->TMPSingle;
	
	VectorField3D* J=new VectorField3D(GS,gridScale,3);
	J->FillValue(0.0);
	float **V=J->V;
	
	if(World->PMF!=NULL)
		for(i=0;i<GS_XYZ;i++)
			UTMP[i]=0.5*(World->IonsQ[ion]*potential[i]+World->PMF[ion][i]);
	else
		for(i=0;i<GS_XYZ;i++)
			UTMP[i]=0.5*World->IonsQ[ion]*potential[i];
	
	fpoh = 4*M_PI*gridScale;
	flcon = 1.6*gridScale*gridScale*1000;
	ftot = flcon/fpoh;
	
	float d,dc,sc;
	double dtmp;
	//Calculate horizontal flew
	int ax;
	for(ax=0;ax<3;ax++)
	{
		switch(ax)
		{
			case 0:
				gridInc=1;
				break;
			case 1:
				gridInc=GS_X;
				break;
			case 2:
				gridInc=GS_XY;
				break;
		}
		ftot = -1.0*World->IonsQ[ion]*flcon/fpoh;
		if(World->D!=NULL)
		{
			for(k=0;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				for(j=0;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=0;i<GS_X-1;i++)
					{
						gridPoint = jgrid+i;
						gridp = gridPoint+gridInc;
						if(NIndexing->GetDiffFloat(ion,gridPoint)>0.0)
						{
							pot = UTMP[gridp]-UTMP[gridPoint];
							d = NIndexing->GetDiffFloat(ion,gridp)>0.0?					0.5*(D[ion][gridPoint]+D[ion][gridp]):0;
							dc = C[ion][gridp]-C[ion][gridPoint];
							sc = C[ion][gridp]+C[ion][gridPoint];
							V[ax][gridPoint]=ftot*d*(dc+pot*sc);
						}
					}
				}
			}
		}
		else
		{
			for(k=1;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				//dtmp=0.0;
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						gridPoint = jgrid+i;
						gridp = gridPoint+gridInc;
						if(NIndexing->GetDiffFloat(ion,gridPoint)>0.0)
						{
							pot = UTMP[gridp]-UTMP[gridPoint];
							d = NIndexing->GetDiffFloat(ion,gridp)>0.0?					0.5*(NIndexing->GetDiffFloat(ion,gridPoint)+NIndexing->GetDiffFloat(ion,gridp)):0;
							dc = C[ion][gridp]-C[ion][gridPoint];
							sc = C[ion][gridp]+C[ion][gridPoint];
							V[ax][gridPoint]=ftot*d*(dc+pot*sc);
						}
					}
				}
			}
		}
	}
	return J;
}
VectorField3D* PoissonNernstPlanckSolver::CalcIinout(VectorField3D* CartI)
{
	int *GS = World->GridSize;
	float gridScale = World->GridScale;
	int GS_X = GS[0];
	int GS_Y = GS[1];
	int GS_Z = GS[2];
	int GS_XY = GS[0]*GS[1];
	int GS_XYZ = GS[2]*GS_XY;
	int ix, iy, iz, GrdPnt;
	VectorField3D* J=new VectorField3D(GS,gridScale,3);
	J->FillValue(0.0);
	for(iz=1;iz<GS_Z-1;iz++)
		for(iy=1;iy<GS_Y-1;iy++)
			for(ix=1;ix<GS_X-1;ix++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		J->V[0][GrdPnt]=0.5*(CartI->V[0][GrdPnt-1]+CartI->V[0][GrdPnt]);
		J->V[1][GrdPnt]=0.5*(CartI->V[1][GrdPnt-GS_X]+CartI->V[1][GrdPnt]);
		J->V[2][GrdPnt]=0.5*(CartI->V[2][GrdPnt-GS_XY]+CartI->V[2][GrdPnt]);
	}
	return J;
}
VectorField3D* PoissonNernstPlanckSolver::CalcAvrI(VectorField3D* I,int iavr)
{
	int *GS = World->GridSize;
	float gridScale = World->GridScale;
	int GS_X = GS[0];
	int GS_Y = GS[1];
	int GS_Z = GS[2];
	int GS_XY = GS[0]*GS[1];
	int GS_XYZ = GS[2]*GS_XY;
	int ix, iy, iz, GrdPnt;
	VectorField3D* J=new VectorField3D(GS,gridScale,3);
	J->FillValue(0.0);
	for(iz=1;iz<GS_Z-1;iz++)
		for(iy=1;iy<GS_Y-1;iy++)
			for(ix=1;ix<GS_X-1;ix++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		J->V[0][GrdPnt]=0.5*(I->V[0][GrdPnt-1]+I->V[0][GrdPnt]);
		J->V[1][GrdPnt]=0.5*(I->V[1][GrdPnt-GS_X]+I->V[1][GrdPnt]);
		J->V[2][GrdPnt]=0.5*(I->V[2][GrdPnt-GS_XY]+I->V[2][GrdPnt]);
	}
	return J;
}
int PoissonNernstPlanckSolver::SolveDouble()
{
	pnpPrint("\nPoisson-Nernst-Plank Solver Double\n");
	int iteration;
	float Tolerance;
	int i,j,k,i1,i2,gridPoint,GrdPnt,ion;
	float gridScale;
	int * gridSize;
	float *PotentialTMP, *Potential;
  //assert(pnpSolverData!=NULL);
#ifdef MPI_PARALLEL
	MPI::Status  status;
#endif
	
	gridScale = World->GridScale;
	gridSize = World->GridSize;
  //PotentialTMP = World->PotentialTMP;
	Potential = World->Potential;
	if(World->PotentialDouble==NULL)
	{
		World->PotentialDouble=new double[World->GS_XYZ];
		for(GrdPnt=0;GrdPnt<World->GS_XYZ;GrdPnt++)
			World->PotentialDouble[GrdPnt]=World->Potential[GrdPnt];
	}
	if(World->CDouble==NULL)
	{
		World->CDouble=new double*[World->NIonsTypes];
		for(ion=0;ion<World->NIonsTypes;ion++)
		{
			World->CDouble[ion]=new double[World->GS_XYZ];
			for(GrdPnt=0;GrdPnt<World->GS_XYZ;GrdPnt++)
				World->CDouble[ion][GrdPnt]=World->C[ion][GrdPnt];
		}
	}
	for(GrdPnt=0;GrdPnt<World->GS_XYZ;GrdPnt++)
		World->PotentialDouble[GrdPnt]=World->Potential[GrdPnt];
	for(ion=0;ion<World->NIonsTypes;ion++)
	{
		for(GrdPnt=0;GrdPnt<World->GS_XYZ;GrdPnt++)
			World->CDouble[ion][GrdPnt]=World->C[ion][GrdPnt];
	}
	
	for(iteration=1;iteration<=MaxIterations;iteration++) 
	{
		int PoissonStatus,NernstPlankStatus;

		assert(Poisson);
		PoissonStatus=Poisson->Solve();
		if(PoissonStatus!=EXIT_SUCCESS)
		{
			return EXIT_FAILURE;
		}
		for(GrdPnt=0;GrdPnt<World->GS_XYZ;GrdPnt++)
			World->PotentialDouble[GrdPnt]=World->Potential[GrdPnt];
		if(PMFWeightMode==1)
		{
			NernstPlank->PMFWeight=double(iteration)/double(MaxIterations);
			pnpPrint0("PMFWeight=%f\n",NernstPlank->PMFWeight);
		}
		NernstPlankStatus=NernstPlank->NernstPlanckSolverDDouble();
    //NernstPlankStatus=NernstPlank->NernstPlanckSolverN2tmp();
		if(Poisson->solver!=Poisson->ArrayDirect)
			Poisson->SetQmobFromConcentrationDouble();
		if(NernstPlankStatus!=EXIT_SUCCESS)
		{
			return EXIT_FAILURE;
		}


		if(World->MyRank==0)
			pnpPrint("<PNPiteration Nit=\"%d\" MaxChange=\"%16.8g\" E=\"%16.8lg\" dE=\"%16.8lg\" rel_dE=\"%16.8lg\"/>\n", iteration, NernstPlank->MaxChange, Poisson->totalEnergy, Poisson->totalChange, Poisson->relativeChange);

		if(iteration%ConvergenceCheck==0)
		{
			for(ion=0;ion<World->NIonsTypes;ion++)
			{
				for(GrdPnt=0;GrdPnt<World->GS_XYZ;GrdPnt++)
					World->C[ion][GrdPnt]=World->CDouble[ion][GrdPnt];
			}
			CombineAndPrintCurrentsDouble(iteration);
		}
	}
	for(GrdPnt=0;GrdPnt<World->GS_XYZ;GrdPnt++)
		World->Potential[GrdPnt]=World->PotentialDouble[GrdPnt];
	for(ion=0;ion<World->NIonsTypes;ion++)
	{
		for(GrdPnt=0;GrdPnt<World->GS_XYZ;GrdPnt++)
			World->C[ion][GrdPnt]=World->CDouble[ion][GrdPnt];
	}
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckSolver::CombineAndPrintCurrents(int iteration)
{
	int i,j,k,i1,i2,gridPoint;
#ifdef MPI_PARALLEL
	MPI::Status  status;
#endif
	int iCurrentDimension;
	iCurrentDimension = currentDimension;
	
	int iLimitCurrentCalcZ[2];
	int diLimitCurrentCalcZ;
	if(bLimitCurrentCalc)
	{
		for(i=0;i<2;i++)
		{
			iLimitCurrentCalcZ[i]=int(World->ConvFloatToGlobIntUnitsZ(LimitCurrentCalcZ[i])+0.5);
			if(iLimitCurrentCalcZ[i]<1)
				iLimitCurrentCalcZ[i]=1;
			if(iLimitCurrentCalcZ[i]>World->GridSizeGlobal[2]-2)
				iLimitCurrentCalcZ[i]=World->GridSizeGlobal[2]-2;
		}
		diLimitCurrentCalcZ=iLimitCurrentCalcZ[1]-iLimitCurrentCalcZ[0]+1;
	}
	
	NernstPlank->nernstPlanckSolverCurrentProfile(iCurrentDimension,positiveCurrentProfile,negativeCurrentProfile);
	float Phi0=0.0,Phi1=0.0;
#ifdef MPI_PARALLEL
	pnpsapp->MyComGroup.Barrier();
	if(World->MyRank==0)
	{
		for(i=1;i<World->NProcs;i++)
		{
			World->GetBorder(&i1,&i2,i);
			pnpsapp->MyComGroup.Recv(positiveCurrentProfile+i1, i2-i1, MPI::DOUBLE, i, 0,status);
			pnpsapp->MyComGroup.Recv(negativeCurrentProfile+i1, i2-i1, MPI::DOUBLE, i, 1,status);
		}
	}
	else
	{
		pnpsapp->MyComGroup.Send(positiveCurrentProfile, World->GridSize[iCurrentDimension]-1, MPI::DOUBLE, 0, 0);
		pnpsapp->MyComGroup.Send(negativeCurrentProfile, World->GridSize[iCurrentDimension]-1, MPI::DOUBLE, 0, 1);
	}
	Phi0=World->Potential[0];
	Phi1=World->Potential[World->GS_XYZ-1];
			
	pnpsapp->MyComGroup.Bcast(&Phi0, 1, MPI::DOUBLE, 0);
	pnpsapp->MyComGroup.Bcast(&Phi1, 1, MPI::DOUBLE, World->NProcs-1);
#else
	DbgPrint0("Phi[0]=%g Phi[%d]=%g MyRank=%d\n",World->Potential[0],World->GS_XYZ-1,World->Potential[World->GS_XYZ-1],World->MyRank);
	Phi0=World->Potential[0];
	Phi1=World->Potential[World->GS_XYZ-1];
#endif
			
	if(World->MyRank==0)
	{
		double currentSum,currentCoordinate;
		double maxCurrent = 0.0;
		double maxCurrentCoordinate = 0.0;
		double currentDev = 0.0;
		double current = 0.0;
		double currentDevPos = 0.0;
		double currentPos = 0.0;
		double currentDevNeg = 0.0;
		double currentNeg = 0.0;

		pnpPrintGroup0("PNPSR-Current =================================================================\n");
		pnpPrintGroup0("PNPSR-Current  iteration %5d\n",iteration);
		pnpPrintGroup0("PNPSR-Current -----------------------------------------------------------------\n");
		pnpPrintGroup0("PNPSR-Current  %8s %16s %16s %16s\n","z, Ang.", "I0, pA", "I1, pA", "Itot, pA");
		pnpPrintGroup0("PNPSR-Current -----------------------------------------------------------------\n");

		//pnpPrint("<PNPCurrent iter=\"%5d\" Format=\"%%6.3f %%16.6f %%16.6f %%16.6f\" Note=\"z Ipos Ineg Itot\">\n",iteration);
		for(k=0;k<World->GridSizeGlobal[iCurrentDimension]-1;k++)
		{
			currentSum = positiveCurrentProfile[k]+negativeCurrentProfile[k];
			pnpPrintGroup0("PNPSR-Current  %8.3f %16.6f %16.6f %16.6f\n", (k-World->GridSizeOriginal[2]/2+0.5)/World->GridScale, positiveCurrentProfile[k], negativeCurrentProfile[k], currentSum);
			//pnpPrint("%6.3f %16.6f %16.6f %16.6f\n",(k-World->GridSizeOriginal[2]/2+0.5)/World->GridScale, positiveCurrentProfile[k], negativeCurrentProfile[k], currentSum);
			currentDev += currentSum*currentSum;
			current += currentSum;
			currentDevPos += positiveCurrentProfile[k]*positiveCurrentProfile[k];
			currentPos += positiveCurrentProfile[k];
			currentDevNeg += negativeCurrentProfile[k]*negativeCurrentProfile[k];
			currentNeg += negativeCurrentProfile[k];
			if(fabs(maxCurrent)<fabs(currentSum)) {
				maxCurrent = currentSum;
				maxCurrentCoordinate = k;
			}
		}
		pnpPrintGroup0("PNPSR-Current -----------------------------------------------------------------\n");

		//pnpPrint("</PNPCurrent>\n");
		currentDev /= World->GridSizeGlobal[iCurrentDimension]-1;
		current /= World->GridSizeGlobal[iCurrentDimension]-1;
		currentDevPos /= World->GridSizeGlobal[iCurrentDimension]-1;
		currentPos /= World->GridSizeGlobal[iCurrentDimension]-1;
		currentDevNeg /= World->GridSizeGlobal[iCurrentDimension]-1;
		currentNeg /= World->GridSizeGlobal[iCurrentDimension]-1;
		ItotErr=currentDev-current*current;
		IposErr=currentDevPos-currentPos*currentPos;
		InegErr=currentDevNeg-currentNeg*currentNeg;
		DbgPrint0("ItotErr=%g IposErr=%g InegErr=%g\n",ItotErr, IposErr, InegErr);
		ItotErr=sqrt(ItotErr);
		IposErr=sqrt(IposErr);
		InegErr=sqrt(InegErr);
		pnpPrintGroup0("PNPSR-Current  maxCurrent = %f pA at z = %f\n",maxCurrent,maxCurrentCoordinate);
		pnpPrintGroup0("PNPSR-Current -----------------------------------------------------------------\n");
		pnpPrintGroup0("PNPSR-Current  Current ( Total ) = %14.7f +- %14.7f pA\n",current, ItotErr);
		pnpPrintGroup0("PNPSR-Current  Current ( Ion0  ) = %14.7f +- %14.7f pA\n",currentPos, IposErr);
		pnpPrintGroup0("PNPSR-Current  Current ( Ion1  ) = %14.7f +- %14.7f pA\n",currentNeg, InegErr);
		pnpPrintGroup0("PNPSR-Current  Applied potential = %14.7f mV\n",(Phi0-Phi1)*CONFAC);
		pnpPrintGroup0("PNPSR-Current -----------------------------------------------------------------\n");
		Itot=current;
		Ipos=currentPos;
		Ineg=currentNeg;
		if(bLimitCurrentCalc)
		{
			currentDev = 0.0;
			current = 0.0;
			currentDevPos = 0.0;
			currentPos = 0.0;
			currentDevNeg = 0.0;
			currentNeg = 0.0;
			for(k=iLimitCurrentCalcZ[0];k<=iLimitCurrentCalcZ[1];k++)
			{
				currentSum = positiveCurrentProfile[k]+negativeCurrentProfile[k];
				currentDev += currentSum*currentSum;
				current += currentSum;
				currentDevPos += positiveCurrentProfile[k]*positiveCurrentProfile[k];
				currentPos += positiveCurrentProfile[k];
				currentDevNeg += negativeCurrentProfile[k]*negativeCurrentProfile[k];
				currentNeg += negativeCurrentProfile[k];
			}
			currentDev /= diLimitCurrentCalcZ;
			current /= diLimitCurrentCalcZ;
			currentDevPos /= diLimitCurrentCalcZ;
			currentPos /= diLimitCurrentCalcZ;
			currentDevNeg /= diLimitCurrentCalcZ;
			currentNeg /= diLimitCurrentCalcZ;
			ItotErr=sqrt(currentDev-current*current);
			IposErr=sqrt(currentDevPos-currentPos*currentPos);
			InegErr=sqrt(currentDevNeg-currentNeg*currentNeg);
			pnpPrintGroup0("PNPSR-Current  Current along Z region: Z = [ %g, %g ]\n", LimitCurrentCalcZ[0], LimitCurrentCalcZ[1]);
			pnpPrintGroup0("PNPSR-Current  CurrentZlimit(Total) = %14.7f +- %14.7f pA\n",current, ItotErr);
			pnpPrintGroup0("PNPSR-Current  CurrentZlimit(Ion0 ) = %14.7f +- %14.7f pA\n",currentPos, IposErr);
			pnpPrintGroup0("PNPSR-Current  CurrentZlimit(Ion1 ) = %14.7f +- %14.7f pA\n",currentNeg, InegErr);
			pnpPrintGroup0("PNPSR-Current  Applied potential    = %14.7f mV\n",(Phi0-Phi1)*CONFAC);
			//pnpPrintGroup0("PNPSR-Current -----------------------------------------------------------------\n");
			//pnpPrint0("<PNPiterCurSummary Zlimit=\"%g %g\" Current=\"%.7e +- %.7e pA\", CurrentPos=\"%.7e +- %.7e pA\" CurrentNeg=\"%.7e +- %.7e pA\" dU=\"%.7e mV\"/>\n",LimitCurrentCalcZ[0],LimitCurrentCalcZ[1],current, ItotErr,currentPos, IposErr,currentNeg, InegErr,(Phi0-Phi1)*CONFAC);
			Itot=current;
			Ipos=currentPos;
			Ineg=currentNeg;
		}
		pnpPrintGroup0("PNPSR-Current =================================================================\n");
	}
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckSolver::CombineAndPrintCurrentsDouble(int iteration)
{
	int i,j,k,i1,i2,gridPoint;
#ifdef MPI_PARALLEL
	MPI::Status  status;
#endif
	int iCurrentDimension;
	iCurrentDimension = currentDimension;
	
	int iLimitCurrentCalcZ[2];
	int diLimitCurrentCalcZ;
	if(bLimitCurrentCalc)
	{
		for(i=0;i<2;i++)
		{
			iLimitCurrentCalcZ[i]=int(World->ConvFloatToGlobIntUnitsZ(LimitCurrentCalcZ[i])+0.5);
			if(iLimitCurrentCalcZ[i]<1)
				iLimitCurrentCalcZ[i]=1;
			if(iLimitCurrentCalcZ[i]>World->GridSizeGlobal[2]-2)
				iLimitCurrentCalcZ[i]=World->GridSizeGlobal[2]-2;
		}
		diLimitCurrentCalcZ=iLimitCurrentCalcZ[1]-iLimitCurrentCalcZ[0]+1;
	}
	
	NernstPlank->nernstPlanckSolverCurrentProfile(iCurrentDimension,positiveCurrentProfile,negativeCurrentProfile);
	float Phi0=0.0,Phi1=0.0;
#ifdef MPI_PARALLEL
	pnpsapp->MyComGroup.Barrier();
	if(World->MyRank==0)
	{
		for(i=1;i<World->NProcs;i++)
		{
			World->GetBorder(&i1,&i2,i);
			pnpsapp->MyComGroup.Recv(positiveCurrentProfile+i1, i2-i1, MPI::DOUBLE, i, 0,status);
			pnpsapp->MyComGroup.Recv(negativeCurrentProfile+i1, i2-i1, MPI::DOUBLE, i, 1,status);
		}
	}
	else
	{
		pnpsapp->MyComGroup.Send(positiveCurrentProfile, World->GridSize[iCurrentDimension]-1, MPI::DOUBLE, 0, 0);
		pnpsapp->MyComGroup.Send(negativeCurrentProfile, World->GridSize[iCurrentDimension]-1, MPI::DOUBLE, 0, 1);
	}
	Phi0=World->Potential[0];
	Phi1=World->Potential[World->GS_XYZ-1];
			
	pnpsapp->MyComGroup.Bcast(&Phi0, 1, MPI::DOUBLE, 0);
	pnpsapp->MyComGroup.Bcast(&Phi1, 1, MPI::DOUBLE, World->NProcs-1);
#else
	DbgPrint0("Phi[0]=%g Phi[%d]=%g MyRank=%d\n",World->Potential[0],World->GS_XYZ-1,World->Potential[World->GS_XYZ-1],World->MyRank);
	Phi0=World->Potential[0];
	Phi1=World->Potential[World->GS_XYZ-1];
#endif
			
	if(World->MyRank==0)
	{
		double currentSum,currentCoordinate;
		double maxCurrent = 0.0;
		double maxCurrentCoordinate = 0.0;
		double currentDev = 0.0;
		double current = 0.0;
		double currentDevPos = 0.0;
		double currentPos = 0.0;
		double currentDevNeg = 0.0;
		double currentNeg = 0.0;
		pnpPrint("<PNPCurrent iter=\"%5d\" Format=\"%%6.3f %%16.6f %%16.6f %%16.6f\" Note=\"z Ipos Ineg Itot\">\n",iteration);
		for(k=0;k<World->GridSizeGlobal[iCurrentDimension]-1;k++)
		{
			currentSum = positiveCurrentProfile[k]+negativeCurrentProfile[k];
			pnpPrint("%6.3f %16.6f %16.6f %16.6f\n",(k-World->GridSizeOriginal[2]/2+0.5)/World->GridScale, positiveCurrentProfile[k], negativeCurrentProfile[k], currentSum);
			currentDev += currentSum*currentSum;
			current += currentSum;
			currentDevPos += positiveCurrentProfile[k]*positiveCurrentProfile[k];
			currentPos += positiveCurrentProfile[k];
			currentDevNeg += negativeCurrentProfile[k]*negativeCurrentProfile[k];
			currentNeg += negativeCurrentProfile[k];
			if(fabs(maxCurrent)<fabs(currentSum)) {
				maxCurrent = currentSum;
				maxCurrentCoordinate = k;
			}
		}
		pnpPrint("</PNPCurrent>\n");
		currentDev /= World->GridSizeGlobal[iCurrentDimension]-1;
		current /= World->GridSizeGlobal[iCurrentDimension]-1;
		currentDevPos /= World->GridSizeGlobal[iCurrentDimension]-1;
		currentPos /= World->GridSizeGlobal[iCurrentDimension]-1;
		currentDevNeg /= World->GridSizeGlobal[iCurrentDimension]-1;
		currentNeg /= World->GridSizeGlobal[iCurrentDimension]-1;
		ItotErr=currentDev-current*current;
		IposErr=currentDevPos-currentPos*currentPos;
		InegErr=currentDevNeg-currentNeg*currentNeg;
		DbgPrint0("ItotErr=%g IposErr=%g InegErr=%g\n",ItotErr, IposErr, InegErr);
		ItotErr=sqrt(ItotErr);
		IposErr=sqrt(IposErr);
		InegErr=sqrt(InegErr);
		pnpPrint0("<pnpSolver: maxCurrent = %f pA at %c = %f\n",maxCurrent,currentDimension,maxCurrentCoordinate);
		pnpPrint0("<PNPiterCurSummary Current=\"%.7e +- %.7e pA\", CurrentPos=\"%.7e +- %.7e pA\" CurrentNeg=\"%.7e +- %.7e pA\" dU=\"%.7e mV\"/>\n",current, ItotErr,currentPos, IposErr,currentNeg, InegErr,(Phi0-Phi1)*CONFAC);
		Itot=current;
		Ipos=currentPos;
		Ineg=currentNeg;
		if(bLimitCurrentCalc)
		{
			currentDev = 0.0;
			current = 0.0;
			currentDevPos = 0.0;
			currentPos = 0.0;
			currentDevNeg = 0.0;
			currentNeg = 0.0;
			for(k=iLimitCurrentCalcZ[0];k<=iLimitCurrentCalcZ[1];k++)
			{
				currentSum = positiveCurrentProfile[k]+negativeCurrentProfile[k];
				currentDev += currentSum*currentSum;
				current += currentSum;
				currentDevPos += positiveCurrentProfile[k]*positiveCurrentProfile[k];
				currentPos += positiveCurrentProfile[k];
				currentDevNeg += negativeCurrentProfile[k]*negativeCurrentProfile[k];
				currentNeg += negativeCurrentProfile[k];
			}
			currentDev /= diLimitCurrentCalcZ;
			current /= diLimitCurrentCalcZ;
			currentDevPos /= diLimitCurrentCalcZ;
			currentPos /= diLimitCurrentCalcZ;
			currentDevNeg /= diLimitCurrentCalcZ;
			currentNeg /= diLimitCurrentCalcZ;
			ItotErr=sqrt(currentDev-current*current);
			IposErr=sqrt(currentDevPos-currentPos*currentPos);
			InegErr=sqrt(currentDevNeg-currentNeg*currentNeg);
			pnpPrint0("<PNPiterCurSummary Zlimit=\"%g %g\" Current=\"%.7e +- %.7e pA\", CurrentPos=\"%.7e +- %.7e pA\" CurrentNeg=\"%.7e +- %.7e pA\" dU=\"%.7e mV\"/>\n",LimitCurrentCalcZ[0],LimitCurrentCalcZ[1],current, ItotErr,currentPos, IposErr,currentNeg, InegErr,(Phi0-Phi1)*CONFAC);
			Itot=current;
			Ipos=currentPos;
			Ineg=currentNeg;
		}
	}
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckSolver::ShowParameters()
{
	cout << "\nParameters of Poisson-Nernst-Plank solver\n";
	cout << "    MaxIterations:.................... "<< MaxIterations <<"\n";
	return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////////////////

PoissonNernstPlanckMultiGridSolver::PoissonNernstPlanckMultiGridSolver()
{
	InitZero();
}


PoissonNernstPlanckMultiGridSolver::~PoissonNernstPlanckMultiGridSolver()
{
	
	Clear();
}
int PoissonNernstPlanckMultiGridSolver::InitZero()
{
	HaObject::SetName("PoissonNernstPlanckMultiGridSolver");
	WorldLev0=NULL;
	PNPLev0=NULL;
	WorldLev1=NULL;
	PNPLev1=NULL;
	
	ItotErr=5550.0;
	IposErr=5550.0;
	InegErr=5550.0;
	positiveCurrentProfile=NULL;
	negativeCurrentProfile=NULL;
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckMultiGridSolver::Clear()
{
	DeleteCArray(positiveCurrentProfile);
	DeleteCArray(negativeCurrentProfile);
	if(PNPLev0!=NULL)
		PNPLev0->NernstPlank->CalcVolume=NULL;
	DeleteObjByPnt(WorldLev0);
	DeleteObjByPnt(PNPLev0);
	DeleteObjByPnt(WorldLev1);
	DeleteObjByPnt(PNPLev1);
	return EXIT_SUCCESS;
}

int PoissonNernstPlanckMultiGridSolver::RunPNPMG(const TiXmlElement* RootElt)
{
	int Status;
	
	const TiXmlElement* CldElt;
	CldElt=RootElt->FirstChildElement();
	do
	{
		if(strcmp("ContWorldLev0",CldElt->Value())==0)
		{
			DeleteObjByPnt(WorldLev0);
			WorldLev0=SingleWorldCalcController::CmdContWorld(CldElt);
			PNP_EXIT_FAIL_NULL(WorldLev0,"Can't initiate ContWorldLev0\n");
		}
		else if(strcmp("ContWorldLev1",CldElt->Value())==0)
		{
			DeleteObjByPnt(WorldLev1);
			WorldLev1=SingleWorldCalcController::CmdContWorld(CldElt);
			PNP_EXIT_FAIL_NULL(WorldLev1,"Can't initiate ContWorldLev1\n");
		}
		else if(strcmp("MapsIODataLev0",CldElt->Value())==0)
		{
			PNP_EXIT_FAIL_NULL(WorldLev0,"WorldLev0 is not initialized\n");
			Status=SingleWorldCalcController::CmdMapIO(CldElt,WorldLev0);
			PNP_EXIT_ON_FAIL_MES(Status,"Can't read maps for ContWorldLev0\n");
		}
		else if(strcmp("MapsIODataLev1",CldElt->Value())==0)
		{
			PNP_EXIT_FAIL_NULL(WorldLev1,"WorldLev1 is not initialized\n");
			Status=SingleWorldCalcController::CmdMapIO(CldElt,WorldLev1);
			PNP_EXIT_ON_FAIL_MES(Status,"Can't read maps for ContWorldLev1\n");
		}
		else if(strcmp("AddPotential",CldElt->Value())==0)
		{
			Status=SingleWorldCalcController::CmdAddPotential(CldElt,WorldLev0);
			PNP_EXIT_ON_FAIL_MES(Status,"Can't Apply Potential on ContWorldLev0\n");
			Status=SingleWorldCalcController::CmdAddPotential(CldElt,WorldLev1);
			PNP_EXIT_ON_FAIL_MES(Status,"Can't Apply Potential on ContWorldLev1\n");
		}
		else if(strcmp("AddDiffPotentialToDiffGroup",CldElt->Value())==0)
		{
			Status=SingleWorldCalcController::CmdAddDiffPotentialToDiffGroup(CldElt,WorldLev0);
			PNP_EXIT_ON_FAIL_MES(Status,"Can't Apply Potential on ContWorldLev0\n");
			Status=SingleWorldCalcController::CmdAddDiffPotentialToDiffGroup(CldElt,WorldLev1);
			PNP_EXIT_ON_FAIL_MES(Status,"Can't Apply Potential on ContWorldLev1\n");
		}
		else if(strcmp("PoissonNernstPlanckSolverLev0",CldElt->Value())==0)
		{
			PNP_EXIT_ON_FAIL_MES(RunPNPonLev0(CldElt),"Can't run PNP on level 0\n");
		}
		else if(strcmp("PoissonNernstPlanckSolverLev1",CldElt->Value())==0)
		{
			PNP_EXIT_ON_FAIL_MES(RunPNPonLev1(CldElt),"Can't run PNP on level 0\n");
		}
		else if(strcmp("PoissonNernstPlanckSolver",CldElt->Value())==0)
		{
			PNP_EXIT_FAIL_NULL(WorldLev0,"WorldLev0 is not initialized\n");
			PNP_EXIT_FAIL_NULL(WorldLev1,"WorldLev1 is not initialized\n");
			if(PNPLev0!=NULL&&PNPLev1!=NULL)
			{
				if(PNPLev0->NernstPlank->CalcVolume==PNPLev0->Poisson->CalcVolume);
					PNPLev0->NernstPlank->CalcVolume=NULL;
			}
			DeleteObjByPnt(PNPLev0);
			DeleteObjByPnt(PNPLev1);
			
			PNP_EXIT_ON_FAIL_MES(ReadPNPnPnNP(CldElt),"Can't read Parameters for PNP run Diffrences\n");
			PNP_EXIT_ON_FAIL_MES(InitSolver(),"Can't init Solver\n");
			PNP_EXIT_ON_FAIL_MES(Solve(),"error during solver run\n");
		}
		else if(strcmp("PoissonSolver",CldElt->Value())==0)
		{
			RunOnlyPoisson(CldElt);
		}
		else if(strcmp("SayHi",CldElt->Value())==0)
		{
			pnpPrint("hi!");
		}
		else
		{
			
		}
		CldElt=CldElt->NextSiblingElement();
	}
	while(CldElt!=NULL);
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckMultiGridSolver::RunPNPonLev0(const TiXmlElement* Elt)
{
	PNP_EXIT_FAIL_NULL(WorldLev0,"WorldLev0 is not initialized\n");
	pnpPrint("will run PNP on level 0\n");
	PoissonNernstPlanckSolver* m_PNPLev0=new PoissonNernstPlanckSolver();
	m_PNPLev0->LoadXML(Elt);
	m_PNPLev0->SetContWorld(WorldLev0);
	m_PNPLev0->InitSolver();
	m_PNPLev0->Solve();
	DeleteObjByPnt(m_PNPLev0)
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckMultiGridSolver::RunPNPonLev1(const TiXmlElement* Elt)
{
	PNP_EXIT_FAIL_NULL(WorldLev1,"WorldLev1 is not initialized\n");
	pnpPrint("will run PNP on level 1\n");
	PoissonNernstPlanckSolver* m_PNPLev1=NULL;
	m_PNPLev1=SingleWorldCalcController::CmdPoissonNernstPlanckSolver(Elt,WorldLev1);
	DeleteObjByPnt(m_PNPLev1)
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckMultiGridSolver::RunOnlyPoisson(const TiXmlElement* Elt)
{
	PoissonSolver *PoissonLev0=new PoissonSolver();
	PoissonLev0->LoadXML(Elt);
	PoissonLev0->SetContWorld(WorldLev0);
	
	PoissonSolver *PoissonLev1=new PoissonSolver();
	PoissonLev1->LoadXML(Elt);
	PoissonLev1->SetContWorld(WorldLev1);
	
	DetermineBoarderBetweenLevels();
	
	//remove points from PNPLev0 which coinside with PNPLev1
	//Poisson-Part
	int i,j,k,GrdPnt;
	
	DeleteCArray(PoissonLev0->CalcVolume);
	if(PoissonLev0->CalcVolume==NULL)
	{
		PoissonLev0->CalcVolume = new bool[WorldLev0->GS_XYZ];
		for(i=1;i<WorldLev0->GS_XYZ;i++)
		{
			PoissonLev0->CalcVolume[i]=false;
		}
		for(k=1;k<WorldLev0->GS_Z-1;k++)
			for(j=1;j<WorldLev0->GS_Y-1;j++)
				for(i=1;i<WorldLev0->GS_X-1;i++)
		{
			GrdPnt = i+j*WorldLev0->GS_X+k*WorldLev0->GS_XY;
			PoissonLev0->CalcVolume[GrdPnt]=true;
		}
		for(k=InternalBoarderMinGridLev0[2];k<=InternalBoarderMaxGridLev0[2];k++)
			for(j=InternalBoarderMinGridLev0[1];j<=InternalBoarderMaxGridLev0[1];j++)
				for(i=InternalBoarderMinGridLev0[0];i<=InternalBoarderMaxGridLev0[0];i++)
		{
			GrdPnt = i+j*WorldLev0->GS_X+k*WorldLev0->GS_XY;
			PoissonLev0->CalcVolume[GrdPnt]=false;
		}
	}
	PoissonLev0->InitSolver();
	PoissonLev1->InitSolver();
	
	int iterP;
	float Tolerance;
	int i1,i2;
	
	int PMaxIterations=PoissonLev0->MaxIterations;
	PoissonLev0->MaxIterations=1;
	PoissonLev1->MaxIterations=1;
	
	//set VectorField3D container for potentials
	float *tmpPotLev0[1];
	tmpPotLev0[0]=WorldLev0->Potential;
	VectorField3D *VFPotLev0=new VectorField3D(WorldLev0->GridSize,WorldLev0->GridScale,1,tmpPotLev0);
	float *tmpPotLev1[1];
	tmpPotLev1[0]=WorldLev1->Potential;
	VectorField3D *VFPotLev1=new VectorField3D(WorldLev1->GridSize,WorldLev1->GridScale,1,tmpPotLev1);
	
	for(iterP=1;iterP<=PMaxIterations;iterP++)
	{
		if(PoissonLev1->PoissonSolverNIB(iterP%PoissonLev1->ConvergenceCheck==0)!=EXIT_SUCCESS)
			return EXIT_FAILURE;
		VFPotLev0->InterpolateInternalBoarderFromExt(VFPotLev1,InternalBoarderMinGridLev0,InternalBoarderMaxGridLev0);
		if(PoissonLev0->PoissonSolverNIB(iterP%PoissonLev1->ConvergenceCheck==0)!=EXIT_SUCCESS)
			return EXIT_FAILURE;
		VFPotLev1->InterpolateBoarderFromExt(VFPotLev0);
	}
		
	delete VFPotLev0;
	delete VFPotLev1;
	delete PoissonLev1;
	delete PoissonLev0;
	
	return EXIT_SUCCESS;
}

int PoissonNernstPlanckMultiGridSolver::ReadPNPnPnNP(const TiXmlElement* Elt)
{
	PNP_EXIT_FAIL_NULL(WorldLev0,"WorldLev0 is not initialized\n");
	PNP_EXIT_FAIL_NULL(WorldLev1,"WorldLev1 is not initialized\n");
	
	PNPLev0=new PoissonNernstPlanckSolver();
	PNPLev0->LoadXML(Elt);
	PNPLev0->SetContWorld(WorldLev0);
	
	PNPLev1=new PoissonNernstPlanckSolver();
	PNPLev1->LoadXML(Elt);
	PNPLev1->SetContWorld(WorldLev1);
	
		//_PoissonNernstPlanckSolver->InitSolver();
		//_PoissonNernstPlanckSolver->Solve();
	if(PNPLev0!=NULL&&PNPLev1!=NULL)
		return EXIT_SUCCESS;
	else
		return EXIT_FAILURE;
}
int PoissonNernstPlanckMultiGridSolver::DetermineBoarderBetweenLevels()
{
	int i;
	int BoarderMinAbsGridLev0[3],BoarderMaxAbsGridLev0[3];
	int BoarderMinAbsGridLev1[3],BoarderMaxAbsGridLev1[3];
	
	float BoarderMinAbsLev0[3],BoarderMaxAbsLev0[3];
	float BoarderMinAbsLev1[3],BoarderMaxAbsLev1[3];
	
	float BrdMaxLev1onLev2AbsGridFloat[3],BrdMinLev1onLev2AbsGridFloat[3];
	int BrdMaxLev1onLev2AbsGrid[3],BrdMinLev1onLev2AbsGrid[3];
	
	int InternalBoarderMinAbsGridLev0[3],InternalBoarderMaxAbsGridLev0[3];
	
	float GridScaleLev0byLev1=WorldLev0->GridScale/WorldLev1->GridScale;
	//int InternalBoarderMinGridLev0[3],InternalBoarderMaxGridLev0[3];
	for(i=0;i<3;i++)
	{
		BoarderMaxAbsGridLev0[i]=WorldLev0->GridSize[i]/2;
		BoarderMinAbsGridLev0[i]=-WorldLev0->GridSize[i]/2;
		BoarderMaxAbsGridLev1[i]=WorldLev1->GridSize[i]/2;
		BoarderMinAbsGridLev1[i]=-WorldLev1->GridSize[i]/2;
		BoarderMaxAbsLev0[i]=BoarderMaxAbsGridLev0[i]/WorldLev0->GridScale;
		BoarderMinAbsLev0[i]=-BoarderMaxAbsLev0[i];
		BoarderMaxAbsLev1[i]=BoarderMaxAbsGridLev1[i]/WorldLev1->GridScale;
		BoarderMinAbsLev1[i]=-BoarderMaxAbsLev1[i];
		
		BrdMaxLev1onLev2AbsGridFloat[i]=BoarderMaxAbsLev1[i]*WorldLev0->GridScale;
		BrdMinLev1onLev2AbsGridFloat[i]=BoarderMinAbsLev1[i]*WorldLev0->GridScale;
		
		BrdMaxLev1onLev2AbsGrid[i]=roundf(BrdMaxLev1onLev2AbsGridFloat[i]);
		BrdMinLev1onLev2AbsGrid[i]=roundf(BrdMinLev1onLev2AbsGridFloat[i]);
		
		InternalBoarderMinAbsGridLev0[i]=roundf(GridScaleLev0byLev1*BoarderMinAbsGridLev1[i])+1;
		InternalBoarderMaxAbsGridLev0[i]=roundf(GridScaleLev0byLev1*BoarderMaxAbsGridLev1[i])-1;
		
		InternalBoarderMinGridLev0[i]=InternalBoarderMinAbsGridLev0[i]+BoarderMaxAbsGridLev0[i];
		InternalBoarderMaxGridLev0[i]=InternalBoarderMaxAbsGridLev0[i]+BoarderMaxAbsGridLev0[i];
	}
	 
	//info
	pnpPrint("\n");
	pnpPrint("WorldLev0: GridSize = [%d, %d, %d], GridScale = %g grid/A\n",WorldLev0->GridSize[0], WorldLev0->GridSize[1], WorldLev0->GridSize[2],WorldLev0->GridScale);
	pnpPrint("           BoxSize = [%g, %g, %g] [A, A, A]\n",WorldLev0->GridSize[0]/WorldLev0->GridScale, WorldLev0->GridSize[1]/WorldLev0->GridScale, WorldLev0->GridSize[2]/WorldLev0->GridScale);
	pnpPrint("           BoarderMinAbsGridLev0 = [%d, %d, %d] [g, g, g]\n",
					 BoarderMinAbsGridLev0[0],BoarderMinAbsGridLev0[1],BoarderMinAbsGridLev0[2]);
	pnpPrint("           BoarderMinAbsLev0 = [%g, %g, %g] [A, A, A]\n",
					 BoarderMinAbsLev0[0],BoarderMinAbsLev0[1],BoarderMinAbsLev0[2]);
	pnpPrint("           BoarderMaxAbsGridLev0 = [%d, %d, %d] [g, g, g]\n",
					 BoarderMaxAbsGridLev0[0],BoarderMaxAbsGridLev0[1],BoarderMaxAbsGridLev0[2]);
	pnpPrint("           BoarderMaxAbsLev0 = [%g, %g, %g] [A, A, A]\n",
					 BoarderMaxAbsLev0[0],BoarderMaxAbsLev0[1],BoarderMaxAbsLev0[2]);
	
	pnpPrint("WorldLev1: GridSize = [%d, %d, %d], GridScale = %g grid/A\n",WorldLev1->GridSize[0], WorldLev1->GridSize[1], WorldLev1->GridSize[2],WorldLev1->GridScale);
	pnpPrint("           BoxSize = [%g, %g, %g] [A, A, A]\n",WorldLev1->GridSize[0]/WorldLev1->GridScale, WorldLev1->GridSize[1]/WorldLev1->GridScale, WorldLev1->GridSize[2]/WorldLev1->GridScale);
	pnpPrint("           BoarderMinAbsGridLev1 = [%d, %d, %d] [g, g, g]\n",
					 BoarderMinAbsGridLev1[0],BoarderMinAbsGridLev1[1],BoarderMinAbsGridLev1[2]);
	pnpPrint("           BoarderMinAbsLev0 = [%g, %g, %g] [A, A, A]\n",
					 BoarderMinAbsLev1[0],BoarderMinAbsLev1[1],BoarderMinAbsLev1[2]);
	pnpPrint("           BoarderMaxAbsGridLev1 = [%d, %d, %d] [g, g, g]\n",
					 BoarderMaxAbsGridLev1[0],BoarderMaxAbsGridLev1[1],BoarderMaxAbsGridLev1[2]);
	pnpPrint("           BoarderMaxAbsLev0 = [%g, %g, %g] [A, A, A]\n",
					 BoarderMaxAbsLev1[0],BoarderMaxAbsLev1[1],BoarderMaxAbsLev1[2]);
	
	pnpPrint("Boarder between levels:\n");
	pnpPrint("           GridScaleLev0byLev1=%g\n",GridScaleLev0byLev1);
	pnpPrint("           BrdMinLev1onLev2AbsGrid = [%d, %d, %d] [g, g, g]\n",
					 BrdMinLev1onLev2AbsGrid[0],BrdMinLev1onLev2AbsGrid[1],BrdMinLev1onLev2AbsGrid[2]);
	pnpPrint("           BrdMinLev1onLev2AbsGridFloat = [%g, %g, %g] [g, g, g]\n",
					 BrdMinLev1onLev2AbsGridFloat[0],BrdMinLev1onLev2AbsGridFloat[1],BrdMinLev1onLev2AbsGridFloat[2]);
	pnpPrint("           BrdMaxLev1onLev2AbsGrid = [%d, %d, %d] [g, g, g]\n",
					 BrdMaxLev1onLev2AbsGrid[0],BrdMaxLev1onLev2AbsGrid[1],BrdMaxLev1onLev2AbsGrid[2]);
	pnpPrint("           BrdMaxLev1onLev2AbsGridFloat = [%g, %g, %g] [g, g, g]\n",
					 BrdMaxLev1onLev2AbsGridFloat[0],BrdMaxLev1onLev2AbsGridFloat[1],BrdMaxLev1onLev2AbsGridFloat[2]);
	pnpPrint("           InternalBoarderMinAbsGridLev0 = [%d, %d, %d] [g, g, g]\n",
					 InternalBoarderMinAbsGridLev0[0],InternalBoarderMinAbsGridLev0[1],InternalBoarderMinAbsGridLev0[2]);
	pnpPrint("           InternalBoarderMinAbsGridLev0 = [%d, %d, %d] [g, g, g]\n",
					 InternalBoarderMaxAbsGridLev0[0],InternalBoarderMaxAbsGridLev0[1],InternalBoarderMaxAbsGridLev0[2]);
	pnpPrint("           InternalBoarderMinGridLev0 = [%d, %d, %d] [g, g, g]\n",
					 InternalBoarderMinGridLev0[0],InternalBoarderMinGridLev0[1],InternalBoarderMinGridLev0[2]);
	pnpPrint("           InternalBoarderMinGridLev0 = [%d, %d, %d] [g, g, g]\n",
					 InternalBoarderMaxGridLev0[0],InternalBoarderMaxGridLev0[1],InternalBoarderMaxGridLev0[2]);
	
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckMultiGridSolver::InitSolver()
{
	DetermineBoarderBetweenLevels();
	
	
	
	//remove points from PNPLev0 which coinside with PNPLev1
	//Poisson-Part
	int i,j,k,GrdPnt;
	
	DeleteCArray(PNPLev0->Poisson->CalcVolume);
	if(PNPLev0->Poisson->CalcVolume==NULL)
	{
		PNPLev0->Poisson->CalcVolume = new bool[WorldLev0->GS_XYZ];
		for(i=1;i<WorldLev0->GS_XYZ;i++)
		{
			PNPLev0->Poisson->CalcVolume[i]=false;
		}
		for(k=1;k<WorldLev0->GS_Z-1;k++)
			for(j=1;j<WorldLev0->GS_Y-1;j++)
				for(i=1;i<WorldLev0->GS_X-1;i++)
		{
			GrdPnt = i+j*WorldLev0->GS_X+k*WorldLev0->GS_XY;
			PNPLev0->Poisson->CalcVolume[GrdPnt]=true;
		}
		for(k=InternalBoarderMinGridLev0[2];k<=InternalBoarderMaxGridLev0[2];k++)
			for(j=InternalBoarderMinGridLev0[1];j<=InternalBoarderMaxGridLev0[1];j++)
				for(i=InternalBoarderMinGridLev0[0];i<=InternalBoarderMaxGridLev0[0];i++)
		{
			GrdPnt = i+j*WorldLev0->GS_X+k*WorldLev0->GS_XY;
			PNPLev0->Poisson->CalcVolume[GrdPnt]=false;
		}
	}
	//Nernst-Planck Part
	PNPLev0->NernstPlank->CalcVolume=PNPLev0->Poisson->CalcVolume;
	PNPLev0->InitSolver();
	PNPLev1->InitSolver();
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckMultiGridSolver::Solve()
{
	//return PNPLev1->Solve();
	pnpPrint("\nPoisson-Nernst-Plank Solver\n");
	int iteration;
	float Tolerance;
	int i,j,k,i1,i2,GrdPnt,IonType;
#ifdef MPI_PARALLEL
	MPI::Status  status;
#endif
	
	int PMaxIterations=PNPLev1->Poisson->MaxIterations;
	int NPMaxIterations=PNPLev1->NernstPlank->MaxIterations;
	PNPLev0->Poisson->MaxIterations=1;
	PNPLev0->NernstPlank->MaxIterations=1;
	PNPLev1->Poisson->MaxIterations=1;
	PNPLev1->NernstPlank->MaxIterations=1;
	int PConvergenceCheck=PNPLev1->Poisson->ConvergenceCheck;
	int NPConvergenceCheck=PNPLev1->NernstPlank->ConvergenceCheck;
	int iterP,iterNP;
	//set VectorField3D container for potentials
	float *tmpPotLev0[1];
	tmpPotLev0[0]=WorldLev0->Potential;
	VectorField3D *VFPotLev0=new VectorField3D(WorldLev0->GridSize,WorldLev0->GridScale,1,tmpPotLev0);
	float *tmpPotLev1[1];
	tmpPotLev1[0]=WorldLev1->Potential;
	VectorField3D *VFPotLev1=new VectorField3D(WorldLev1->GridSize,WorldLev1->GridScale,1,tmpPotLev1);
	//set VectorField3D container for concentration
	VectorField3D *VFConcLev0=new VectorField3D(WorldLev0->GridSize,WorldLev0->GridScale,WorldLev0->NIonsTypes,WorldLev0->C);
	VectorField3D *VFConcLev1=new VectorField3D(WorldLev1->GridSize,WorldLev1->GridScale,WorldLev1->NIonsTypes,WorldLev1->C);
	
	//float fpohLev0= 4*M_PI*WorldLev1->GridScale;
	float coefLev0 = COANGS / (WorldLev0->GridScale*WorldLev0->GridScale*WorldLev0->GridScale);
	//float fpohLev1= 4*M_PI*WorldLev1->GridScale;
	float coefLev1 = COANGS / (WorldLev1->GridScale*WorldLev1->GridScale*WorldLev1->GridScale);
	float NPMult=coefLev0/coefLev1;
	
	DbgPrint0("NPMult=%5g\n",NPMult);
	
	for(iteration=1;iteration<=PNPLev1->MaxIterations;iteration++) 
	{
		
		//Poisson in MG
		for(iterP=1;iterP<=PMaxIterations;iterP++)
		{
			if(PNPLev1->Poisson->PoissonSolverNIB(iterP==PMaxIterations||iterP%PConvergenceCheck==0)!=EXIT_SUCCESS)
				return EXIT_FAILURE;
			VFPotLev0->InterpolateInternalBoarderFromExt(VFPotLev1,InternalBoarderMinGridLev0,InternalBoarderMaxGridLev0);
			if(PNPLev0->Poisson->PoissonSolverNIB(iterP==PMaxIterations||iterP%PConvergenceCheck==0)!=EXIT_SUCCESS)
				return EXIT_FAILURE;
			VFPotLev1->InterpolateBoarderFromExt(VFPotLev0);
		}
		for(IonType=0;IonType<WorldLev0->NIonsTypes;IonType++)
		{
			PNPLev0->NernstPlank->NernstPlanckSolverDpre(IonType);
			PNPLev1->NernstPlank->NernstPlanckSolverDpre(IonType);
			for(iterNP=1;iterNP<=NPMaxIterations;iterNP++)
			{
				if(PNPLev1->NernstPlank->NernstPlanckSolverDiter(IonType,iterNP==NPMaxIterations||iterNP%NPConvergenceCheck==0)!=EXIT_SUCCESS)
					return EXIT_FAILURE;
				VFConcLev0->InterpolateInternalBoarderFromExtWithMultipl(VFConcLev1,InternalBoarderMinGridLev0,InternalBoarderMaxGridLev0,NPMult);
				
				if(PNPLev0->NernstPlank->NernstPlanckSolverDiter(IonType,iterNP==NPMaxIterations||iterNP%NPConvergenceCheck==0)!=EXIT_SUCCESS)
					return EXIT_FAILURE;
				VFConcLev1->InterpolateBoarderFromExtWithMultipl(VFConcLev0,1.0/NPMult);
			}
			PNPLev0->NernstPlank->NernstPlanckSolverDpost(IonType);
			PNPLev1->NernstPlank->NernstPlanckSolverDpost(IonType);
		}
// 		for(i=0;i<WorldLev0->GS_XYZ;i++)
// 		{
// 			if(PNPLev0->Poisson->CalcVolume[i])
// 			WorldLev0->Potential[i]=2.0;
// 		}
		/*for(k=InternalBoarderMinGridLev0[2];k<=InternalBoarderMaxGridLev0[2];k++)
			for(j=InternalBoarderMinGridLev0[1];j<=InternalBoarderMaxGridLev0[1];j++)
				for(i=InternalBoarderMinGridLev0[0];i<=InternalBoarderMaxGridLev0[0];i++)
		{
			GrdPnt = i+j*WorldLev0->GS_X+k*WorldLev0->GS_XY;
			WorldLev0->Potential[GrdPnt]=0.0;
		}*/
// 		NernstPlankStatus=0;
// 		if(PNPLev1->NernstPlank->Solve()==EXIT_SUCCESS)
// 			NernstPlankStatus++;
// 		if(PNPLev0->NernstPlank->Solve()==EXIT_SUCCESS)
// 			NernstPlankStatus++;
// 		if(NernstPlankStatus!=2)
// 		{
// 			return EXIT_FAILURE;
// 		}
		
		if(PNPLev1->Poisson->solver!=PNPLev1->Poisson->ArrayDirect)
			PNPLev1->Poisson->SetQmobFromConcentration();
		
		if(PNPLev0->Poisson->solver!=PNPLev0->Poisson->ArrayDirect)
			PNPLev0->Poisson->SetQmobFromConcentration();


		if(WorldLev1->MyRank==0)
			pnpPrint("<PNPiterationLev1 Nit=\"%d\" MaxChange=\"%16.8g\" E=\"%16.8lg\" dE=\"%16.8lg\" rel_dE=\"%16.8lg\"/>\n", iteration, PNPLev1->NernstPlank->MaxChange, PNPLev1->Poisson->totalEnergy, PNPLev1->Poisson->totalChange, PNPLev1->Poisson->relativeChange);
		if(WorldLev1->MyRank==0)
			pnpPrint("<PNPiterationLev0 Nit=\"%d\" MaxChange=\"%16.8g\" E=\"%16.8lg\" dE=\"%16.8lg\" rel_dE=\"%16.8lg\"/>\n", iteration, PNPLev0->NernstPlank->MaxChange, PNPLev0->Poisson->totalEnergy, PNPLev0->Poisson->totalChange, PNPLev0->Poisson->relativeChange);
		if(iteration%PNPLev1->ConvergenceCheck==0)
		{
			//PNPLev0->CombineAndPrintCurrents(iteration);
			PNPLev1->CombineAndPrintCurrents(iteration);
		}
	}
	PNPLev0->Poisson->MaxIterations=PMaxIterations;
	PNPLev0->NernstPlank->MaxIterations=NPMaxIterations;
	PNPLev1->Poisson->MaxIterations=PMaxIterations;
	PNPLev1->NernstPlank->MaxIterations=NPMaxIterations;
	
	delete VFPotLev0;
	delete VFPotLev1;
	delete VFConcLev0;
	delete VFConcLev1;
	return EXIT_SUCCESS;
}
