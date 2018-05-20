//
// C++ Implementation: poissonsolver
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

#include "poissonsolver.h"
#include "tinyxml.h"
#include "pnpdebug.h"
#include "contworld.h"
#include "math.h"
#include "pnpconstants.h"
#include "mapio.h"

#include "pmfcalculation.h"
#include "pnputil.h"
#define FMUL(X,Y) (X*Y)
#define FADD(X,Y) (X+Y)
#define FMAF(X,Y,Z) (X*Y+Z)

#include <stdlib.h>

PoissonSolver::PoissonSolver()
 : GenericSolver()
{
	InitZero();
}


PoissonSolver::~PoissonSolver()
{
	Clear();
}
int PoissonSolver::InitZero()
{
	HaObject::SetName("PoissonSolver");
	SolverStr.push_back("Auto");
	SolverStr.push_back("NodeIndexBased");
	SolverStr.push_back("ArrayDirect");
	SolverStr.push_back("PNPC");
	
	World=NULL;
	
	ConvergenceCheck=20;
	QmobMod=0;
	MaxIterations=300;
	Convergence=0.0;
	Relaxation=1.9;
	solver=0;
	verbose=true;
	QmobMod=2;
	MinIterations=0;

	 //Initial values of variables:
	totalChange=0;
	relativeChange=0;
	totalEnergy=0;
	totalEnergyInd=0;
	verbose=true;
	
	dielectricZSSUM=NULL;
	dielectricZDBSUM=NULL;
	
	dielectricXS = NULL;
	dielectricYS = NULL;
	dielectricZS = NULL;
	dielectricZSSUM = NULL;
	dielectricXmS = NULL;
	dielectricYmS = NULL;
	dielectricZmS = NULL;
	QstS = NULL;
	dielectricXDB = NULL;
	dielectricYDB = NULL;
	dielectricZDB = NULL;
	dielectricZDBSUM = NULL;
	dielectricXmDB = NULL;
	dielectricYmDB = NULL;
	dielectricZmDB = NULL;
	dielectricCh = NULL;
	
	IndexNoSingular=NULL;
	IndexDielBoarder=NULL;
	IndexCharge=NULL;
	IndexSingular=NULL;
	//ChargeSum=NULL;
	Qst=NULL;
	PhiCharge=NULL;
	PhiSingular=NULL;
	QmobNum[0]=0;
	QmobNum[1]=0;
	QmobNum[2]=0;
	IndexQmob=NULL;
	dielectricChMob=NULL;
	Qmob=NULL;
	
	QmobDielBoarderNum[0]=0;
	QmobDielBoarderNum[1]=0;
	QmobDielBoarderNum[2]=0;
	IndexQmobDielBoarder=NULL;
	QmobDielBoarder=NULL;
	dielectricXQmobDB=NULL;
	dielectricYQmobDB=NULL;
	dielectricZQmobDB=NULL;
	dielectricXmQmobDB=NULL;
	dielectricYmQmobDB=NULL;
	dielectricZmQmobDB=NULL;
	dielectricZQmobDBSUM=NULL;
	
	QmobDielBoarderQstNum[0]=0;
	QmobDielBoarderQstNum[1]=0;
	QmobDielBoarderQstNum[2]=0;
	IndexQmobDielBoarderQst=NULL;
	QmobDielBoarderQst=NULL;
	QstQmobDielBoarderQst=NULL;
	dielectricXQmobDBQst=NULL;
	dielectricYQmobDBQst=NULL;
	dielectricZQmobDBQst=NULL;
	dielectricXmQmobDBQst=NULL;
	dielectricYmQmobDBQst=NULL;
	dielectricZmQmobDBQst=NULL;
	dielectricZQmobDBSUMQst=NULL;
	
	//QmobFlag=NULL;
	
	ChargeSum=NULL;
	
	ChargesEnergy=NULL;
	potential=NULL;
	
	WayToCalcSystemEnergy=0;
	ConvFacMaxHistory=1;

	CalcVolume=NULL;
	return EXIT_SUCCESS;
}
int PoissonSolver::Clear()
{
//	DeleteCArray(IonsQ);
//	DeleteCVecArray(C,NIonsTypes);
//	DeleteObjByPnt(NIndexing);
	DeleteCArray(IndexQmob);
	DeleteCArray(Qmob);
	//DeleteCArray(QmobFlag);
	DeleteCArray(dielectricChMob);
	DeleteCArray(dielectricXS);
	DeleteCArray(dielectricYS);
	DeleteCArray(dielectricZS);
	DeleteCArray(dielectricZSSUM);
	DeleteCArray(dielectricXmS);
	DeleteCArray(dielectricYmS);
	DeleteCArray(dielectricZmS);
	DeleteCArray(QstS);
	DeleteCArray(dielectricXDB);
	DeleteCArray(dielectricYDB);
	DeleteCArray(dielectricZDB);
	DeleteCArray(dielectricZDBSUM);
	DeleteCArray(dielectricXmDB);
	DeleteCArray(dielectricYmDB);
	DeleteCArray(dielectricZmDB);
	DeleteCArray(dielectricCh);
	DeleteCArray(IndexNoSingular);
	DeleteCArray(IndexDielBoarder);
	DeleteCArray(IndexCharge);
	DeleteCArray(IndexSingular);
	DeleteCArray(Qst);
	
	DeleteCArray(IndexQmobDielBoarder);
	DeleteCArray(QmobDielBoarder);
	DeleteCArray(dielectricXQmobDB);
	DeleteCArray(dielectricYQmobDB);
	DeleteCArray(dielectricZQmobDB);
	DeleteCArray(dielectricXmQmobDB);
	DeleteCArray(dielectricYmQmobDB);
	DeleteCArray(dielectricZmQmobDB);
	DeleteCArray(dielectricZQmobDBSUM);
	
	DeleteCArray(IndexQmobDielBoarderQst);
	DeleteCArray(QmobDielBoarderQst);
	DeleteCArray(QstQmobDielBoarderQst);
	DeleteCArray(dielectricXQmobDBQst);
	DeleteCArray(dielectricYQmobDBQst);
	DeleteCArray(dielectricZQmobDBQst);
	DeleteCArray(dielectricXmQmobDBQst);
	DeleteCArray(dielectricYmQmobDBQst);
	DeleteCArray(dielectricZmQmobDBQst);
	DeleteCArray(dielectricZQmobDBSUMQst);
	
	DeleteCArray(ChargesEnergy);
	DeleteCArray(PhiCharge);
	DeleteCArray(PhiSingular);
	
	DeleteCArray(CalcVolume);
#	ifdef PNPDOUBLE
	DeleteCArray(potential);
#	endif
	return EXIT_SUCCESS;
}
int PoissonSolver::SaveXML(TiXmlElement* Elt, HaContext* p_ctxt )
{
	return EXIT_SUCCESS;
}
int PoissonSolver::LoadXML(const TiXmlElement* Elt, HaContext* p_ctxt )
{
	if(strncmp(HaObject::GetCStrName(),Elt->Value(),13))
	{
		fprintf(stderr,"ERROR: Wrong XML Element %s, expecting %s\n",Elt->Value(),HaObject::GetCStrName());
		return EXIT_FAILURE;
	}
	Clear();
	int i,gridPoint;
	//Read Primary Parameters
	Elt->GetIntAttribute("MaxIterations",&MaxIterations);
	Elt->GetFloatAttribute("Convergence",&Convergence);
	Elt->GetFloatAttribute("Relaxation",&Relaxation);  
	if(Elt->GetIntAttribute("ConvergenceCheck",&ConvergenceCheck)!=EXIT_SUCCESS)ConvergenceCheck=20;
	if(Elt->GetStdStrIndex("Solver",&solver,SolverStr)!=EXIT_SUCCESS)solver=0;
	if(Elt->GetBoolAttribute("Verbose",&verbose)!=EXIT_SUCCESS)verbose=true;
	if(Elt->GetIntAttribute("QmobMod",&QmobMod)!=EXIT_SUCCESS)QmobMod=2;
	if(Elt->GetIntAttribute("MinIterations",&MinIterations)!=EXIT_SUCCESS)MinIterations=0;
	
	if(Elt->GetIntAttribute("ConvFacMaxHistory",&ConvFacMaxHistory)!=EXIT_SUCCESS)ConvFacMaxHistory=1;
	ShowParameters();
	
	SetRelaxation(Relaxation);
	
	if(dielectricZSSUM!=NULL)ShowProperties();
	
	return EXIT_SUCCESS;
}
int PoissonSolver::LoadParamFromPyDict(PyObject *dict)
{
	Clear();
	int i,gridPoint;
	//Read Primary Parameters
	haPyDict_GetItemAsInt(dict,"MaxIterations",&MaxIterations);
	haPyDict_GetItemAsFloat(dict,"Relaxation",&Relaxation);
	haPyDict_GetItemAsFloat(dict,"Convergence",&Convergence);
	haPyDict_GetItemAsFloat(dict,"Tolerance",&Convergence);
	haPyDict_GetItemAsInt(dict,"ConvergenceCheck",&ConvergenceCheck);
	haPyDict_GetItemAsInt(dict,"Solver",&solver);
	haPyDict_GetItemAsBool(dict,"Verbose",&verbose);

	haPyDict_GetItemAsInt(dict,"QmobMod",&QmobMod);
	haPyDict_GetItemAsInt(dict,"MinIterations",&MinIterations);
	haPyDict_GetItemAsInt(dict,"ConvFacMaxHistory",&ConvFacMaxHistory);

	ShowParameters();
	
	SetRelaxation(Relaxation);
	
	if(dielectricZSSUM!=NULL)ShowProperties();
	
	return EXIT_SUCCESS;
}
int PoissonSolver::SetRelaxation(float _Relaxation)
{
	int i,gridPoint;
	Relaxation=_Relaxation;
	//Scale Parameters to Internal Units
	//Calculate Derivative Parameters
	if(dielectricZSSUM!=NULL)for(i=0;i<SingularNum[2];i++){
		gridPoint=IndexSingular[i];
		dielectricZSSUM[i]=dielectricXS[i]+dielectricYS[i]+dielectricZS[i]+dielectricXmS[i]+dielectricYmS[i]+dielectricZmS[i];
		dielectricZSSUM[i]=Relaxation/dielectricZSSUM[i];
	}
	if(dielectricZDBSUM!=NULL)for(i=0;i<DielBoarderNum[2];i++){
		gridPoint=IndexDielBoarder[i];
		dielectricZDBSUM[i]=dielectricXDB[i]+dielectricYDB[i]+dielectricZDB[i]+dielectricXmDB[i]+dielectricYmDB[i]+dielectricZmDB[i];
		dielectricZDBSUM[i]=Relaxation/dielectricZDBSUM[i];
	}
	if(dielectricZQmobDBSUM!=NULL)for(i=0;i<QmobDielBoarderNum[2];i++){
		gridPoint=IndexQmobDielBoarder[i];
		dielectricZQmobDBSUM[i]=dielectricXQmobDB[i]+dielectricYQmobDB[i]+dielectricZQmobDB[i]+dielectricXmQmobDB[i]+dielectricYmQmobDB[i]+dielectricZmQmobDB[i];
		dielectricZQmobDBSUM[i]=Relaxation/dielectricZQmobDBSUM[i];
	}
	om2 = Relaxation;
	om1 = 1.0-om2;
	om2d6 = om2/6.0;
	
	if(dielectricZSSUM!=NULL)ShowProperties();
	return EXIT_SUCCESS;
}
int PoissonSolver::SetContWorld(ContWorld* _world)
{
	int i;
	World=_world;
	GridScale = World->GridScale;
	GS_X = World->GridSize[0];
	GS_Y = World->GridSize[1];
	GS_Z = World->GridSize[2];
	GS_XY = GS_X*GS_Y;
	GS_XYZ = GS_XY*GS_Z;
	if(World->Potential == NULL){
		if(!(World->Potential = new float[GS_XYZ])){
			fprintf(stderr,"ERROR 104: No memory available\n");
			exit(104);
		}
		for(i=0;i<GS_XYZ;i++)World->Potential[i]=0.0;
	}
#	ifndef PNPDOUBLE
	potential = World->Potential;
#	else
	
#	endif
	return EXIT_SUCCESS;
}
int PoissonSolver::ShowParameters()
{
	DbgPrint2("CPoisson::ShowParameters\n");
	
	pnpPrintGroup0("\nParameters of Poisson solver set up\n");
	pnpPrintGroup0("    MaxIterations:.................... %d\n", MaxIterations);
	pnpPrintGroup0("    Convergence:...................... %.8g kT\n", Convergence);
	pnpPrintGroup0("    Relaxation:....................... %.5g\n", Relaxation);
	
	return EXIT_SUCCESS;
}
int PoissonSolver::ShowProperties()
{
	DbgPrint2("CPoisson::ShowProperties\n");
	return EXIT_SUCCESS;
}
int PoissonSolver::InitSolver()
{
  //Solver{Auto=0,NodeIndexBased=1,ArrayDirect=2,PNPC=3};
  //Solver solver;
	int status;
	bool bGuessNumberOfIteration=false;
	if(Relaxation<0.0||MaxIterations<0)
	{
		bGuessNumberOfIteration=true;
		Relaxation=1.0;
	}
	if(solver==Auto)
	{
		if(World->NIndexing==NULL)status=InitSolverAD();
		else status=InitSolverNIB();
	}
	else if(solver==NodeIndexBased)
		status=InitSolverNIB();
	else if(solver==ArrayDirect)
		status=InitSolverAD();
	else if(solver==PNPC)
		status=InitSolverW();
	
	int itmp;
	if(bGuessNumberOfIteration)
		itmp=GuessNumberOfIteration();
	if(MaxIterations<0)
	{
		MaxIterations=itmp*abs(MaxIterations);
		if(MaxIterations<60)
			MaxIterations=60;
	}
	return status;
}
int PoissonSolver::InitSolverW()
{
	DbgPrint0("CPoisson::InitSolverW()");
	long count = 0;
	int gridSizeX;
	int gridSizeY;
	int gridSizeZ;
	long gridSizeXY;
	long gridPoint;
	int i,j,k;
	float * dielectric[3];
	float * staticCharge;
	long jgrid,kgrid;

//   assert(output!=NULL);
//   assert(output->dielectricMap[0]!=NULL);
//   assert(output->dielectricMap[1]!=NULL);
//   assert(output->dielectricMap[2]!=NULL);
//   assert(output->potentialMap!=NULL);
//   assert(output->staticChargeMap!=NULL);
	gridSizeX = World->GridSize[0];
	gridSizeY = World->GridSize[1];
	gridSizeZ = World->GridSize[2];
	gridSizeXY = gridSizeX*gridSizeY;
	
	poissonSolverData = new PoissonSolverDataW();
	assert(poissonSolverData!=NULL);
	
	poissonSolverData->maxIterations = MaxIterations;
	poissonSolverData->convergence = Convergence;
	poissonSolverData->relaxation = Relaxation;
	
	World->CheckArrays("P",true);
	
	if(World->C!=NULL)
		poissonSolverData->hasDynamicCharges = 1;
	else
		poissonSolverData->hasDynamicCharges = 0;
	
	
	for(i=0;i<3;i++)
	{
		dielectric[i] = World->Epsilon[i];
	}
	staticCharge = World->Qstat;

	for(k=1;k<gridSizeZ-1;k++) {
		kgrid = k*gridSizeXY;
		for(j=1;j<gridSizeY-1;j++) {
			jgrid = kgrid+j*gridSizeX;
			for(i=1;i<gridSizeX-1;i++) {
				gridPoint = jgrid+i;
				if(dielectric[0][gridPoint]!=dielectric[0][gridPoint-1] || dielectric[0][gridPoint]!=dielectric[1][gridPoint] || dielectric[0][gridPoint]!=dielectric[1][gridPoint-gridSizeX] || dielectric[0][gridPoint]!=dielectric[2][gridPoint] || dielectric[0][gridPoint]!=dielectric[2][gridPoint-gridSizeXY] || (poissonSolverData->hasDynamicCharges==0 && staticCharge[gridPoint]!=0))
					count++;
			}
		}
	}

	poissonSolverData->borderPoints = new long[count+1];
	assert(poissonSolverData->borderPoints!=NULL);
	poissonSolverData->om2InverseDielectricSum = new float [count];
	assert(poissonSolverData->om2InverseDielectricSum!=NULL);
	if(!poissonSolverData->hasDynamicCharges) {
		poissonSolverData->typeOfBorderPoint = new short[count];
		assert(poissonSolverData->typeOfBorderPoint!=NULL);
	}
	else
		poissonSolverData->typeOfBorderPoint = NULL;

	count = 0;

	for(k=1;k<gridSizeZ-1;k++) {
		kgrid = k*gridSizeXY;
		for(j=1;j<gridSizeY-1;j++) {
			jgrid = kgrid+j*gridSizeX;
			for(i=1;i<gridSizeX-1;i++) {
				gridPoint = jgrid+i;
				if(dielectric[0][gridPoint]!=dielectric[0][gridPoint-1] || dielectric[0][gridPoint]!=dielectric[1][gridPoint] || dielectric[0][gridPoint]!=dielectric[1][gridPoint-gridSizeX] || dielectric[0][gridPoint]!=dielectric[2][gridPoint] || dielectric[0][gridPoint]!=dielectric[2][gridPoint-gridSizeXY]) {
					if(poissonSolverData->hasDynamicCharges==0) {
						if(staticCharge[gridPoint]!=0)
							poissonSolverData->typeOfBorderPoint[count] = POISSON_SOLVER_CHARGED_UNCOMMON_DIELECTRIC_POINT;
						else 
							poissonSolverData->typeOfBorderPoint[count] = POISSON_SOLVER_UNCOMMON_DIELECTRIC_POINT;
					}
					poissonSolverData->om2InverseDielectricSum[count] = Relaxation/(dielectric[0][gridPoint]+dielectric[0][gridPoint-1]+dielectric[1][gridPoint]+dielectric[1][gridPoint-gridSizeX]+dielectric[2][gridPoint]+dielectric[2][gridPoint-gridSizeXY]);
					poissonSolverData->borderPoints[count] = gridPoint;
					count++;
				}
				else if(poissonSolverData->hasDynamicCharges==0 && staticCharge[gridPoint]!=0) {
					poissonSolverData->typeOfBorderPoint[count] = POISSON_SOLVER_CHARGED_POINT;
					poissonSolverData->borderPoints[count] = gridPoint;
					count++;
				}
			}
		}
	}
	poissonSolverData->borderPoints[count] = gridSizeXY*gridSizeZ; 
	return EXIT_SUCCESS;
}
int PoissonSolver::InitSolverNIB()
{
	DbgPrint2("CPoisson::InitSolverNIB\n");
	if(World==NULL)return NULL;
	
	//if(World->C!=NULL)SetQmobForPNP();
	//temp vars
	int IType;
	int i,j,k,kgrid,jgrid,BlackOrWhite;
	int GrdPnt;
	
	NodeIndex* NIndex=World->NIndexing->NIndex;
  
	//If volume there run calculation is not set, make a fake one
//	Clear();
	
	//count b/w regions
	NoSingularNum[0]=0;NoSingularNum[1]=0;NoSingularNum[2]=0;
	SingularNum[0]=0;SingularNum[1]=0;SingularNum[2]=0;
	DielBoarderNum[0]=0;DielBoarderNum[1]=0;DielBoarderNum[2]=0;
	ChargeNum[0]=0;ChargeNum[1]=0;ChargeNum[2]=0;
	
	
	//correction for spesific type of calculation, P(no additional charge) P(NP)(dynamic charge)
	unsigned int specChargeMask=NodeIndexing::ChargeMask;
	//!@todo not forget to change when will do pnp
//   if(World->D!=NULL)//i.e. P(NP)
//   {
//     specChargeMask=ChargeMask|DiffMask;
//   }
	unsigned int ChargeDielBoarderMask=specChargeMask|NodeIndexing::DielBoarderMask;
	unsigned int BlackAndWhiteMask=NodeIndexing::BlackAndWhiteMask;
	unsigned int DielBoarderMask=NodeIndexing::DielBoarderMask;
	DbgPrint0("DielBoarderMask=%X\n",DielBoarderMask);
	DbgPrint0("ChargeDielBoarderMask=%X\n",ChargeDielBoarderMask);
	DbgPrint0("BlackAndWhiteMask=%X\n",BlackAndWhiteMask);
	DbgPrint0("specChargeMask=%X\n",specChargeMask);
	
	int QmobTot=QmobNum[2]+QmobDielBoarderNum[2]+QmobDielBoarderQstNum[2];
	bool bQmobHere;
	bool bCalcVolume;
	bCalcVolume=true;
	if(QmobTot==0)
	{
		for(k=1;k<GS_Z-1;k++)
		{
			kgrid = k*GS_XY;
			for(j=1;j<GS_Y-1;j++)
			{
				jgrid = kgrid+j*GS_X;
				for(i=1;i<GS_X-1;i++)
				{
					GrdPnt = jgrid+i;
					if(CalcVolume!=NULL)
						bCalcVolume=CalcVolume[GrdPnt];
					if(bCalcVolume)
					{
						if((NIndex[GrdPnt]&ChargeDielBoarderMask)==ChargeDielBoarderMask)
						{
							SingularNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
							SingularNum[2]++;
						}
						else if(NIndex[GrdPnt]&specChargeMask)
						{
							ChargeNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
							ChargeNum[2]++;
						}
						else if(NIndex[GrdPnt]&DielBoarderMask)
						{
							DielBoarderNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
							DielBoarderNum[2]++;
						}
						else
						{
							NoSingularNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
							NoSingularNum[2]++;
						}
					}
				}
			}
		}
	}
	else
	{
		int CountQmobOnDielBoarder=0;
		int CountQmobOnQst=0;
		int CountQmobOnDielBoarderAndQst=0;
		for(k=1;k<GS_Z-1;k++)
		{
			kgrid = k*GS_XY;
			for(j=1;j<GS_Y-1;j++)
			{
				jgrid = kgrid+j*GS_X;
				for(i=1;i<GS_X-1;i++)
				{
					GrdPnt = jgrid+i;
					if(CalcVolume!=NULL)
						bCalcVolume=CalcVolume[GrdPnt];
					if(bCalcVolume)
					{
						bQmobHere=false;
						for(IType=0;IType<World->NIonsTypes;IType++)
							if(World->C[IType][GrdPnt]>0.0)bQmobHere=true;
						if(bQmobHere==true)
						{
							if((NIndex[GrdPnt]&ChargeDielBoarderMask)==ChargeDielBoarderMask)
								CountQmobOnDielBoarderAndQst++;
							else
							{
								if(NIndex[GrdPnt]&DielBoarderMask)CountQmobOnDielBoarder++;
								if(NIndex[GrdPnt]&specChargeMask)CountQmobOnQst++;
							}
						}
						else if((NIndex[GrdPnt]&ChargeDielBoarderMask)==ChargeDielBoarderMask)
						{
							
							SingularNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
							SingularNum[2]++;
						}
						else if(NIndex[GrdPnt]&specChargeMask)
						{
							ChargeNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
							ChargeNum[2]++;
						}
						else if(NIndex[GrdPnt]&DielBoarderMask)
						{
							//DbgPrint0("NIndex[%d]=%X\n",GrdPnt,NIndex[GrdPnt]);
							unsigned int t1=NIndex[GrdPnt];
							DielBoarderNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
							DielBoarderNum[2]++;
						}
						else
						{
							NoSingularNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
							NoSingularNum[2]++;
						}
					}
				}
			}
		}
		DbgPrint0("CountQmobOnDielBoarder=%d CountQmobOnQst=%d CountQmobOnDielBoarderAndQst=%d\n",CountQmobOnDielBoarder, CountQmobOnQst, CountQmobOnDielBoarderAndQst);
		if(CountQmobOnQst>0)
		{
			fprintf(stderr,"ERROR\nERROR situation then Qmob On Qst is not implemented yet\n");
		}
	}
	DbgPrint1("InitSolver:Total: MyrankWorld->MyRank=%d icharge=% d    \niDielBoarder=%d iSingular=%d iNoSingular=%d QmobNum=%d\n",
						World->MyRank,ChargeNum[2],DielBoarderNum[2],SingularNum[2],NoSingularNum[2],QmobNum[2]);
	DbgPrint0("TotalNodes: %d\n",ChargeNum[2]+DielBoarderNum[2]+SingularNum[2]+NoSingularNum[2]+QmobNum[2]);
	DbgPrint1("InitSolver:Blach cell:: MyrankWorld->MyRank=%d icharge=% d    \niDielBoarder=%d iSingular=%d iNoSingular=%d QmobNum=%d\n",
						World->MyRank,ChargeNum[1],DielBoarderNum[1],SingularNum[1],NoSingularNum[1],QmobNum[2]);
  
	if(!(IndexNoSingular = new int[NoSingularNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(IndexDielBoarder = new int[DielBoarderNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(IndexCharge = new int[ChargeNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(Qst = new float[ChargeNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(dielectricCh = new float[ChargeNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(IndexSingular = new int[SingularNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(QstS = new float[SingularNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(PhiSingular = new float[SingularNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(PhiCharge = new float[ChargeNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	for(i=0;i<ChargeNum[2];i++)
		PhiCharge[i]=0.0f;
	for(i=0;i<SingularNum[2];i++)
		PhiSingular[i]=0.0f;
	
  //fill indexes of nodes types
	int iCharge=0,iDielBoarder=0,iSingular=0,iNoSingular=0;
	int iCharge2=ChargeNum[1],iDielBoarder2=DielBoarderNum[1],iSingular2=SingularNum[1],iNoSingular2=NoSingularNum[1];
	int QCount=0;
	float *Q=World->NIndexing->Q;
	float *Eps=World->NIndexing->Eps;
	double q=0.0f;
	for(i=0;i<World->NIndexing->QNum;i++)
	{
		q+=Q[i];
	}
	DbgPrint0("q=%f QCount=%d GridScale=%f\n",(float)q/4/M_PI/World->GridScale,World->NIndexing->QNum,World->GridScale);
	q=0.0f;
	if(QmobTot==0)
	{
		for(k=1;k<GS_Z-1;k++)
		{
			kgrid = k*GS_XY;
			for(j=1;j<GS_Y-1;j++)
			{
				jgrid = kgrid+j*GS_X;
				for(i=1;i<GS_X-1;i++)
				{
					GrdPnt = jgrid+i;
					if(CalcVolume!=NULL)
						bCalcVolume=CalcVolume[GrdPnt];
					if(bCalcVolume)
					{
						if((NIndex[GrdPnt]&ChargeDielBoarderMask)==ChargeDielBoarderMask)
						{
							if(NIndex[GrdPnt]&BlackAndWhiteMask)
							{
								IndexSingular[iSingular]=GrdPnt;
								QstS[iSingular]=Q[QCount];
								q+=Q[QCount];
								iSingular++;
							}
							else
							{
								IndexSingular[iSingular2]=GrdPnt;
								QstS[iSingular2]=Q[QCount];
								q+=Q[QCount];
								iSingular2++;
							}
						}
						else if(NIndex[GrdPnt]&specChargeMask)
						{
							if(NIndex[GrdPnt]&BlackAndWhiteMask)
							{
								IndexCharge[iCharge]=GrdPnt;
								dielectricCh[iCharge]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
								Qst[iCharge]=Q[QCount]/dielectricCh[iCharge];
								q+=Q[QCount];
								iCharge++;
							}
							else
							{
								IndexCharge[iCharge2]=GrdPnt;
								dielectricCh[iCharge2]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
								Qst[iCharge2]=Q[QCount]/dielectricCh[iCharge2];
								q+=Q[QCount];
								iCharge2++;
							}
						}
						else if(NIndex[GrdPnt]&DielBoarderMask)
						{
							if(NIndex[GrdPnt]&BlackAndWhiteMask)
							{
								IndexDielBoarder[iDielBoarder]=GrdPnt;
								iDielBoarder++;
							}
							else
							{
								IndexDielBoarder[iDielBoarder2]=GrdPnt;
								iDielBoarder2++;
							}
						}
						else
						{
							if(NIndex[GrdPnt]&BlackAndWhiteMask)
							{
								IndexNoSingular[iNoSingular]=GrdPnt;
								iNoSingular++;
							}
							else
							{
								IndexNoSingular[iNoSingular2]=GrdPnt;
								iNoSingular2++;
							}
						}
						if(NIndex[GrdPnt]&specChargeMask)
							QCount++;
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
			for(j=1;j<GS_Y-1;j++)
			{
				jgrid = kgrid+j*GS_X;
				for(i=1;i<GS_X-1;i++)
				{
					GrdPnt = jgrid+i;
					if(CalcVolume!=NULL)
						bCalcVolume=CalcVolume[GrdPnt];
					if(bCalcVolume)
					{
						bQmobHere=false;
						for(IType=0;IType<World->NIonsTypes;IType++)
							if(World->C[IType][GrdPnt]>0.0)bQmobHere=true;
						if(bQmobHere==true)
						{
							
						}
						else if((NIndex[GrdPnt]&ChargeDielBoarderMask)==ChargeDielBoarderMask)
						{
							if(NIndex[GrdPnt]&BlackAndWhiteMask)
							{
								IndexSingular[iSingular]=GrdPnt;
								QstS[iSingular]=Q[QCount];
								q+=Q[QCount];
								iSingular++;
							}
							else
							{
								IndexSingular[iSingular2]=GrdPnt;
								QstS[iSingular2]=Q[QCount];
								q+=Q[QCount];
								iSingular2++;
							}
						}
						else if(NIndex[GrdPnt]&specChargeMask)
						{
							if(NIndex[GrdPnt]&BlackAndWhiteMask)
							{
								IndexCharge[iCharge]=GrdPnt;
								dielectricCh[iCharge]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
								Qst[iCharge]=Q[QCount]/dielectricCh[iCharge];
								q+=Q[QCount];
								iCharge++;
							}
							else
							{
								IndexCharge[iCharge2]=GrdPnt;
								dielectricCh[iCharge2]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
								Qst[iCharge2]=Q[QCount]/dielectricCh[iCharge2];
								q+=Q[QCount];
								iCharge2++;
							}
						}
						else if(NIndex[GrdPnt]&DielBoarderMask)
						{
							if(NIndex[GrdPnt]&BlackAndWhiteMask)
							{
								IndexDielBoarder[iDielBoarder]=GrdPnt;
								iDielBoarder++;
							}
							else
							{
								IndexDielBoarder[iDielBoarder2]=GrdPnt;
								iDielBoarder2++;
							}
						}
						else
						{
							if(NIndex[GrdPnt]&BlackAndWhiteMask)
							{
								IndexNoSingular[iNoSingular]=GrdPnt;
								iNoSingular++;
							}
							else
							{
								IndexNoSingular[iNoSingular2]=GrdPnt;
								iNoSingular2++;
							}
						}
						if(NIndex[GrdPnt]&specChargeMask)
							QCount++;
					}
				}
			}
		}
	}
	DbgPrint0("q=%f QCount=%d",(float)q/4/M_PI/World->GridScale,QCount);
  //Allocate dielectric helping array
	dielectricXDB = new float[DielBoarderNum[2]];
	dielectricYDB = new float[DielBoarderNum[2]];
	dielectricZDB = new float[DielBoarderNum[2]];
	dielectricXmDB = new float[DielBoarderNum[2]];
	dielectricYmDB = new float[DielBoarderNum[2]];
	dielectricZmDB = new float[DielBoarderNum[2]];
	dielectricZDBSUM = new float[DielBoarderNum[2]];
	dielectricXS = new float[SingularNum[2]];
	dielectricYS = new float[SingularNum[2]];
	dielectricZS = new float[SingularNum[2]];
	dielectricXmS = new float[SingularNum[2]];
	dielectricYmS = new float[SingularNum[2]];
	dielectricZmS = new float[SingularNum[2]];
	dielectricZSSUM = new float[SingularNum[2]];
  

  
	if(!dielectricXDB||!dielectricYDB||!dielectricZDB||!dielectricZDBSUM||!dielectricXS||!dielectricYS||!dielectricZS||!dielectricZSSUM){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!dielectricXmDB||!dielectricYmDB||!dielectricZmDB||!dielectricXmS||!dielectricYmS||!dielectricZmS){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
  
  //float *Epsilon[3];
  //for(i=0;i<3;i++)Epsilon[i] = World->Epsilon[i];
  
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		dielectricXS[i]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
		dielectricYS[i]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
		dielectricZS[i]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
		dielectricXmS[i]=Eps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
		dielectricYmS[i]=Eps[(NIndex[GrdPnt-GS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
		dielectricZmS[i]=Eps[(NIndex[GrdPnt-GS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
		dielectricZSSUM[i]=dielectricXS[i]+dielectricYS[i]+dielectricZS[i]+dielectricXmS[i]+dielectricYmS[i]+dielectricZmS[i];
		dielectricZSSUM[i]=Relaxation/dielectricZSSUM[i];
	}
	for(i=0;i<DielBoarderNum[2];i++)
	{
		GrdPnt=IndexDielBoarder[i];
		dielectricXDB[i]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
		dielectricYDB[i]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
		dielectricZDB[i]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
		dielectricXmDB[i]=Eps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
		dielectricYmDB[i]=Eps[(NIndex[GrdPnt-GS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
		dielectricZmDB[i]=Eps[(NIndex[GrdPnt-GS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
		dielectricZDBSUM[i]=dielectricXDB[i]+dielectricYDB[i]+dielectricZDB[i]+dielectricXmDB[i]+dielectricYmDB[i]+dielectricZmDB[i];
		dielectricZDBSUM[i]=Relaxation/dielectricZDBSUM[i];
	}
  
  /*VectorField3D* Vf=new VectorField3D(World->GridSize,World->GridScale,1);
	float *v=Vf->V[0];
	j=0;
	for(i=0;i<World->GridSizeXYZ;i++)v[i]=0.0;
	for(i=SingularNum[0];i<SingularNum[2];i++){
	GrdPnt=IndexSingular[i];
	v[GrdPnt]+=1.0;
	j++;
}
	for(i=ChargeNum[0];i<ChargeNum[2];i++){
	GrdPnt=IndexCharge[i];
	v[GrdPnt]+=1.0;
	j++;
}
  
	for(i=DielBoarderNum[0];i<DielBoarderNum[2];i++){
	GrdPnt=IndexDielBoarder[i];
	v[GrdPnt]+=1.0;
	j++;
       //potential[GrdPnt]+=100.0;
}
	for(i=NoSingularNum[0];i<NoSingularNum[2];i++){
	GrdPnt=IndexNoSingular[i];
	v[GrdPnt]+=1.0;
	j++;
}
	for(i=QmobNum[0];i<QmobNum[2];i++){
	GrdPnt=IndexQmob[i];
	v[GrdPnt]+=1.0;
	j++;
}
	DbgPrint0("j=%d\n QmobNum[2]=%d",j,QmobNum[2]);
	Vf->WriteToFile("map.bin");
	delete Vf;*/
	return EXIT_SUCCESS;
}
int PoissonSolver::InitSolverAD()
{
      /// @todo implement me
  //std::istrstream iStr(ProcedureCommand);
  //char tmp[MAP_IO_STRING_LENGTH];
  
  //Read parameters
	DbgPrint2("Renew Singularities List\n");

	int gridSizeX;
	int gridSizeY;
	int gridSizeZ;
	int gridSizeXY,gridSizeXYZ;
	int gridPoint;
	int i,j,k;
	float * dielectric[3];
	float * staticCharge;
	int jgrid,kgrid;
	int iCharge,iDielBoarder,iSingular,iNoSingular;
	int BlackOrWhite;
	float *positiveCharge;
	float *negativeCharge;
  
	int *typeOfBorderPoint=NULL;
	
	World->CheckArrays("P",true);
	
	if(World->C==NULL)
	{
		positiveCharge = NULL;
		negativeCharge = NULL;
	}
	else
	{
		positiveCharge = World->C[0];
		negativeCharge = World->C[1];
	}
	gridSizeX = World->GridSize[0];
	gridSizeY = World->GridSize[1];
	gridSizeZ = World->GridSize[2];
	gridSizeXY = gridSizeX*gridSizeY;
	gridSizeXYZ = gridSizeZ*gridSizeXY;

	if(dielectricXS!=NULL){delete [] dielectricXS;dielectricXS=NULL;}
	if(dielectricYS!=NULL){delete [] dielectricYS;dielectricYS=NULL;}
	if(dielectricZS!=NULL){delete [] dielectricZS;dielectricZS=NULL;}
	if(dielectricZSSUM!=NULL){delete [] dielectricZSSUM;dielectricZSSUM=NULL;}
	if(dielectricXmS!=NULL){delete [] dielectricXmS;dielectricXmS=NULL;}
	if(dielectricYmS!=NULL){delete [] dielectricYmS;dielectricYmS=NULL;}
	if(dielectricZmS!=NULL){delete [] dielectricZmS;dielectricZmS=NULL;}
	if(dielectricXDB!=NULL){delete [] dielectricXDB;dielectricXDB=NULL;}
	if(dielectricYDB!=NULL){delete [] dielectricYDB;dielectricYDB=NULL;}
	if(dielectricZDB!=NULL){delete [] dielectricZDB;dielectricZDB=NULL;}
	if(dielectricZDBSUM!=NULL){delete [] dielectricZDBSUM;dielectricZDBSUM=NULL;}
	if(dielectricXmDB!=NULL){delete [] dielectricXmDB;dielectricXmDB=NULL;}
	if(dielectricYmDB!=NULL){delete [] dielectricYmDB;dielectricYmDB=NULL;}
	if(dielectricZmDB!=NULL){delete [] dielectricZmDB;dielectricZmDB=NULL;}
	if(dielectricCh!=NULL){delete [] dielectricCh;dielectricCh=NULL;}
	if(IndexNoSingular!=NULL){delete [] IndexNoSingular;IndexNoSingular=NULL;}
	if(IndexDielBoarder!=NULL){delete [] IndexDielBoarder;IndexDielBoarder=NULL;}
	if(IndexCharge!=NULL){delete [] IndexCharge;IndexCharge=NULL;}
	if(IndexSingular!=NULL){delete [] IndexSingular;IndexSingular=NULL;}
	if(ChargeSum!=NULL){delete [] ChargeSum;ChargeSum=NULL;}
  

  //Allocate Potential
	if(World->Potential == NULL){
		if(!(World->Potential = new float[gridSizeXYZ])){
			fprintf(stderr,"ERROR 104: No memory available\n");
			exit(104);
		}
		for(j=0;j<GS_XYZ;j++)World->Potential[j]=0.0;
	}
  

	for(i=0;i<3;i++)
		dielectric[i] = World->Epsilon[i];
  
	if(!(typeOfBorderPoint = new int[GS_XYZ])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}


	if(!(ChargeSum = new float[GS_XYZ])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	} 
	staticCharge=ChargeSum;
  
	if(QmobMod==0||World->D==NULL)
	{
		for(i=0;i<GS_XYZ;i++)
			if(World->Qstat[i]!=0)staticCharge[i]=1.0;
	}
	else {
		for(i=0;i<GS_XYZ;i++)
			if(World->Qstat[i]!=0||World->D[0][i]>0||World->D[1][i]>0)staticCharge[i]=1.0;
	}
	for(k=0;k<GS_XYZ;k++){
    //dielectricSum[k] = 6*dielectric[0][k];
		typeOfBorderPoint[k] = Boarder;
	}

  
	NoSingularNum[0]=0;NoSingularNum[1]=0;NoSingularNum[2]=0;
	SingularNum[0]=0;SingularNum[1]=0;SingularNum[2]=0;
	DielBoarderNum[0]=0;DielBoarderNum[1]=0;DielBoarderNum[2]=0;
	ChargeNum[0]=0;ChargeNum[1]=0;ChargeNum[2]=0;
  
  //int *Temp=new int[gridSizeXYZ];
	for(k=1;k<gridSizeZ-1;k++) {
		kgrid = k*gridSizeXY;
		for(j=1;j<gridSizeY-1;j++) {
			jgrid = kgrid+j*gridSizeX;
			for(i=1;i<gridSizeX-1;i++) {
				gridPoint = jgrid+i;
				BlackOrWhite=k+j+i+World->startBlackAndWhite;
				typeOfBorderPoint[gridPoint] = NoSingular;
        
        //if(BlackOrWhite%2==0)Temp[gridPoint]=0;
        //else Temp[gridPoint]=gridPoint;
				if(dielectric[0][gridPoint]!=dielectric[0][gridPoint-1] ||
							 dielectric[0][gridPoint]!=dielectric[1][gridPoint]||
							 dielectric[0][gridPoint]!=dielectric[1][gridPoint-gridSizeX] ||
							 dielectric[0][gridPoint]!=dielectric[2][gridPoint] ||
							 dielectric[0][gridPoint]!=dielectric[2][gridPoint-gridSizeXY])
				{
          //dielectricSum[gridPoint] = dielectric[0][gridPoint] + dielectric[0][gridPoint-1] + dielectric[1][gridPoint] + dielectric[1][gridPoint-gridSizeX] + dielectric[2][gridPoint] + dielectric[2][gridPoint-gridSizeXY];
					if(staticCharge[gridPoint]!=0){
						typeOfBorderPoint[gridPoint] = ChargeAndDielBoarder;
						if(BlackOrWhite%2==0)SingularNum[1]++;
						SingularNum[2]++;
					}
					else{
						typeOfBorderPoint[gridPoint] = DielBoarder;
						if(BlackOrWhite%2==0)DielBoarderNum[1]++;
						DielBoarderNum[2]++;
					}
				}
				else if(staticCharge[gridPoint]!=0) {
					typeOfBorderPoint[gridPoint] = Charge;
					if(BlackOrWhite%2==0)ChargeNum[1]++;
					ChargeNum[2]++;
				}
				if(typeOfBorderPoint[gridPoint] == NoSingular){
					if(BlackOrWhite%2==0)NoSingularNum[1]++;
					NoSingularNum[2]++;
				}
			}
		}
	}
  //if(World->MyRank==0)World->writeMapFloat("tmpde0.gz",World->Epsilon[1],gridSizeXYZ,1);
  
	DbgPrint1("Total: MyrankWorld->MyRank=%d icharge=% d    \niDielBoarder=%d iSingular=%d iNoSingular=%d\n",
						World->MyRank,ChargeNum[2],DielBoarderNum[2],SingularNum[2],NoSingularNum[2]);
	DbgPrint1("Blach cell:: MyrankWorld->MyRank=%d icharge=% d    \niDielBoarder=%d iSingular=%d iNoSingular=%d\n",
						World->MyRank,ChargeNum[1],DielBoarderNum[1],SingularNum[1],NoSingularNum[1]);
    
	if(!(IndexNoSingular = new int[NoSingularNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(IndexDielBoarder = new int[DielBoarderNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(IndexCharge = new int[ChargeNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(IndexSingular = new int[SingularNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(PhiSingular = new float[SingularNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(PhiCharge = new float[ChargeNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	for(i=0;i<ChargeNum[2];i++)
		PhiCharge[i]=0.0f;
	for(i=0;i<SingularNum[2];i++)
		PhiSingular[i]=0.0f;
	
	iNoSingular=-1;
	iSingular=-1; 
	iCharge=-1;
	iDielBoarder=-1; 
    //for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint=gridPoint+2){
	for(k=1;k<gridSizeZ-1;k++) {
		kgrid = k*gridSizeXY;
		for(j=1;j<gridSizeY-1;j++) {
			jgrid = kgrid+j*gridSizeX;
			for(i=1;i<gridSizeX-1;i++) {
				gridPoint = jgrid+i;
				BlackOrWhite=k+j+i+World->startBlackAndWhite;
				if(BlackOrWhite%2==0){
					switch(typeOfBorderPoint[gridPoint]){
						case NoSingular:
							iNoSingular++;
							IndexNoSingular[iNoSingular]=gridPoint;
							break;
						case Charge:
							iCharge++;
							IndexCharge[iCharge]=gridPoint;
							break;
						case DielBoarder:
							iDielBoarder++;
							IndexDielBoarder[iDielBoarder]=gridPoint;
							break;
						case ChargeAndDielBoarder:
							iSingular++;
							IndexSingular[iSingular]=gridPoint;
							break;
						case Boarder:
							break;   
					}
				}
			}
		}
	}
    //for(gridPoint=1;gridPoint<gridSizeXYZ;gridPoint=gridPoint+2){
	for(k=1;k<gridSizeZ-1;k++) {
		kgrid = k*gridSizeXY;
		for(j=1;j<gridSizeY-1;j++) {
			jgrid = kgrid+j*gridSizeX;
			for(i=1;i<gridSizeX-1;i++) {
				gridPoint = jgrid+i;
				BlackOrWhite=k+j+i+World->startBlackAndWhite;
				if(BlackOrWhite%2==1){
					switch(typeOfBorderPoint[gridPoint]){
						case NoSingular:
							iNoSingular++;
							IndexNoSingular[iNoSingular]=gridPoint;
							break;
						case Charge:
							iCharge++;
							IndexCharge[iCharge]=gridPoint;
							break;
						case DielBoarder:
							iDielBoarder++;
							IndexDielBoarder[iDielBoarder]=gridPoint;
							break;
						case ChargeAndDielBoarder:
							iSingular++;
							IndexSingular[iSingular]=gridPoint;
							break;
						case Boarder:
							break;   
					}
				}
			}
		}
	}
  //Free memory. SAVE it for children!
	if(typeOfBorderPoint!=NULL)
	{
		delete [] typeOfBorderPoint;
		typeOfBorderPoint=NULL;
	}
  
	//
  
	dielectricXDB = new float[DielBoarderNum[2]];
	dielectricYDB = new float[DielBoarderNum[2]];
	dielectricZDB = new float[DielBoarderNum[2]];
	dielectricXmDB = new float[DielBoarderNum[2]];
	dielectricYmDB = new float[DielBoarderNum[2]];
	dielectricZmDB = new float[DielBoarderNum[2]];
	dielectricZDBSUM = new float[DielBoarderNum[2]];
	dielectricXS = new float[SingularNum[2]];
	dielectricYS = new float[SingularNum[2]];
	dielectricZS = new float[SingularNum[2]];
	dielectricXmS = new float[SingularNum[2]];
	dielectricYmS = new float[SingularNum[2]];
	dielectricZmS = new float[SingularNum[2]];
	dielectricZSSUM = new float[SingularNum[2]];
  
	if(!(dielectricCh = new float[ChargeNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
  
	if(!dielectricXDB||!dielectricYDB||!dielectricZDB||!dielectricZDBSUM||!dielectricXS||!dielectricYS||!dielectricZS||!dielectricZSSUM){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!dielectricXmDB||!dielectricYmDB||!dielectricZmDB||!dielectricXmS||!dielectricYmS||!dielectricZmS){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
  
	for(i=0;i<SingularNum[2];i++){
		gridPoint=IndexSingular[i];
		dielectricXS[i]=dielectric[0][gridPoint];
		dielectricYS[i]=dielectric[1][gridPoint];
		dielectricZS[i]=dielectric[2][gridPoint];
		dielectricXmS[i]=dielectric[0][gridPoint-1];
		dielectricYmS[i]=dielectric[1][gridPoint-gridSizeX];
		dielectricZmS[i]=dielectric[2][gridPoint-gridSizeXY];
		dielectricZSSUM[i]=dielectricXS[i]+dielectricYS[i]+dielectricZS[i]+dielectricXmS[i]+dielectricYmS[i]+dielectricZmS[i];
		dielectricZSSUM[i]=Relaxation/dielectricZSSUM[i];
	}
	for(i=0;i<ChargeNum[2];i++){
		gridPoint=IndexCharge[i];
    //denominator[gridPoint]=temp1;
		dielectricCh[i]=dielectric[0][gridPoint];
    //staticCharge[gridPoint]/=dielectricCh[i];
	}
	for(i=0;i<DielBoarderNum[2];i++){
		gridPoint=IndexDielBoarder[i];
		dielectricXDB[i]=dielectric[0][gridPoint];
		dielectricYDB[i]=dielectric[1][gridPoint];
		dielectricZDB[i]=dielectric[2][gridPoint];
		dielectricXmDB[i]=dielectric[0][gridPoint-1];
		dielectricYmDB[i]=dielectric[1][gridPoint-gridSizeX];
		dielectricZmDB[i]=dielectric[2][gridPoint-gridSizeXY];
		dielectricZDBSUM[i]=dielectricXDB[i]+dielectricYDB[i]+dielectricZDB[i]+dielectricXmDB[i]+dielectricYmDB[i]+dielectricZmDB[i];
		dielectricZDBSUM[i]=Relaxation/dielectricZDBSUM[i];
	}
//   if(MemoryCarefullUsage)
//   {
//     delete [] World->Epsilon[0];
//     World->Epsilon[0]=NULL;
//     delete [] World->Epsilon[1];
//     World->Epsilon[1]=NULL;
//     delete [] World->Epsilon[2];
//     World->Epsilon[2]=NULL;
//   }
	return EXIT_SUCCESS;
}
int PoissonSolver::SetQmobFromConcentration()
{
	int i, IType;
	int GrdPnt;
	float **C=World->C;
	double q;
	
	for(i=QmobNum[0];i<QmobNum[2];i++)
	{
		GrdPnt=IndexQmob[i];
		Qmob[i]=(C[0][GrdPnt]-C[1][GrdPnt])/dielectricChMob[i];
	}
	for(i=QmobDielBoarderNum[0];i<QmobDielBoarderNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarder[i];
		QmobDielBoarder[i]=(C[0][GrdPnt]-C[1][GrdPnt]);
	}
	for(i=QmobDielBoarderQstNum[0];i<QmobDielBoarderQstNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarderQst[i];
		QmobDielBoarderQst[i]=(C[0][GrdPnt]-C[1][GrdPnt]);
	}
	q=0.0;
	for(i=0;i<QmobNum[2];i++)
	{
		q+=Qmob[i]*dielectricChMob[i];
	}
	for(i=0;i<QmobDielBoarderNum[2];i++)
	{
		q+=QmobDielBoarder[i];
	}
	for(i=0;i<QmobDielBoarderQstNum[2];i++)
	{
		q+=QmobDielBoarderQst[i];
	}
	DbgPrint0("Mobile Charges induced: q=%f QmobNum=%d\n",(float)q/4/M_PI/World->GridScale);
	return EXIT_SUCCESS;
}
int PoissonSolver::SetQmobFromConcentrationDouble()
{
	int i, IType;
	int GrdPnt;
	double **C=World->CDouble;
	double q;
	
	for(i=QmobNum[0];i<QmobNum[2];i++)
	{
		GrdPnt=IndexQmob[i];
		q=(C[0][GrdPnt]-C[1][GrdPnt])/dielectricChMob[i];
		Qmob[i]=q;
	}
	for(i=QmobDielBoarderNum[0];i<QmobDielBoarderNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarder[i];
		q=(C[0][GrdPnt]-C[1][GrdPnt]);
		QmobDielBoarder[i]=q;
	}
	for(i=QmobDielBoarderQstNum[0];i<QmobDielBoarderQstNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarderQst[i];
		q=(C[0][GrdPnt]-C[1][GrdPnt]);
		QmobDielBoarderQst[i]=q;
	}
	q=0.0;
	for(i=0;i<QmobNum[2];i++)
	{
		q+=Qmob[i]*dielectricChMob[i];
	}
	for(i=0;i<QmobDielBoarderNum[2];i++)
	{
		q+=QmobDielBoarder[i];
	}
	for(i=0;i<QmobDielBoarderQstNum[2];i++)
	{
		q+=QmobDielBoarderQst[i];
	}
	DbgPrint0("Mobile Charges induced: q=%f QmobNum=%d\n",(float)q/4/M_PI/World->GridScale);
	return EXIT_SUCCESS;
}
int PoissonSolver::SetQmobForPNP()
{
	int i,j,k, i1, IType;
	int ix, iy, iz;
	int GrdPnt;
	int jgrid,kgrid;
	int GS_X = World->GridSize[0];
	int GS_Y = World->GridSize[1];
	int GS_Z = World->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	int itmp1,itmp2;
	int iQmob,iQmob2;
	int iQmobDielBoarder,iQmobDielBoarder2;
	int iQmobDielBoarderQst,iQmobDielBoarderQst2;
	
	NodeIndexing* NIndexing=World->NIndexing;
	NodeIndex* NIndex=NIndexing->NIndex;
	
	NodeIndex DiffBoarderMask=NodeIndexing::DiffIon0BoarderMask;
	NodeIndex BlackAndWhiteMask=NodeIndexing::BlackAndWhiteMask;
	NodeIndex DielBoarderMask=NodeIndexing::DielBoarderMask;
	
	NodeIndex specChargeMask=NodeIndexing::ChargeMask;
	NodeIndex ChargeDielBoarderMask=specChargeMask|NodeIndexing::DielBoarderMask;
	
	PNP_EXIT_FAIL_NULL(World->C,"Concentration Maps do not exist\n");
	for(IType=0;IType<World->NIonsTypes;IType++)
	{
		PNP_EXIT_FAIL_NULL(World->C[IType],"One of Concentration Maps does not exist\n");
	}
	PNP_EXIT_FAIL_NULL(World->NIndexing,"World->NIndexing does not exist\n");
	
	if(World->D==NULL)
	{
		for(IType=0;IType<World->NIonsTypes;IType++)
		{
			for(i=0;i<GS_XYZ;i++)
			{
					if(NIndexing->GetDiffFloat(IType,i)==0.0)World->C[IType][i]=0.0f;
			}
		}
	}
	else
	{
		for(IType=0;IType<World->NIonsTypes;IType++)
		{
			for(i=0;i<GS_XYZ;i++)
			{
				if(World->D[IType][i]==0.0)World->C[IType][i]=0.0f;
			}
		}
	}
	QmobNum[0]=0;
	QmobNum[1]=0;
	QmobNum[2]=0;
	QmobDielBoarderNum[0]=0;
	QmobDielBoarderNum[1]=0;
	QmobDielBoarderNum[2]=0;
	int countChargeDielBoarder=0;
	
	bool bCalcVolume=true;
	bool bQmobHere;
	//!@todo dising for 2 ions type
	for(k=1;k<GS_Z-1;k++)
	{
		kgrid = k*GS_XY;
		for(j=1;j<GS_Y-1;j++)
		{
			jgrid = kgrid+j*GS_X;
			for(i=1;i<GS_X-1;i++)
			{
				GrdPnt = jgrid+i;
				if(CalcVolume!=NULL)
					bCalcVolume=CalcVolume[GrdPnt];
				
				bQmobHere=false;
				for(IType=0;IType<World->NIonsTypes;IType++)
					if(World->C[IType][GrdPnt]>0.0)bQmobHere=true;
				
				if(bCalcVolume&&bQmobHere)
				{
					if((NIndex[GrdPnt]&ChargeDielBoarderMask)==ChargeDielBoarderMask)
					{
						QmobDielBoarderQstNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
						QmobDielBoarderQstNum[2]++;
					}
					else if(NIndex[GrdPnt]&DielBoarderMask)
					{
						QmobDielBoarderNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
						QmobDielBoarderNum[2]++;
					}
					else
					{
						QmobNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
						QmobNum[2]++;
					}
					
				}
			}
		}
	}
	if(countChargeDielBoarder>0)
		pnpError("countChargeDielBoarder=%d>0 this situation is not imlemented\n",countChargeDielBoarder);
	if(QmobNum[2]>0)
	{
		if(IndexQmob==NULL)IndexQmob=new int[QmobNum[2]];
		if(Qmob==NULL)Qmob=new float[QmobNum[2]];
		if(dielectricChMob==NULL)dielectricChMob=new float[QmobNum[2]];
	}
	if(QmobDielBoarderNum[2]>0)
	{
		if(IndexQmobDielBoarder==NULL)IndexQmobDielBoarder=new int[QmobDielBoarderNum[2]];
		if(QmobDielBoarder==NULL)QmobDielBoarder=new float[QmobDielBoarderNum[2]];
		dielectricXQmobDB = new float[QmobDielBoarderNum[2]];
		dielectricYQmobDB = new float[QmobDielBoarderNum[2]];
		dielectricZQmobDB = new float[QmobDielBoarderNum[2]];
		dielectricXmQmobDB = new float[QmobDielBoarderNum[2]];
		dielectricYmQmobDB = new float[QmobDielBoarderNum[2]];
		dielectricZmQmobDB = new float[QmobDielBoarderNum[2]];
		dielectricZQmobDBSUM = new float[QmobDielBoarderNum[2]];
	}
	if(QmobDielBoarderQstNum[2]>0)
	{
		if(IndexQmobDielBoarderQst==NULL)IndexQmobDielBoarderQst=new int[QmobDielBoarderQstNum[2]];
		if(QmobDielBoarderQst==NULL)QmobDielBoarderQst=new float[QmobDielBoarderQstNum[2]];
		if(QstQmobDielBoarderQst==NULL)QstQmobDielBoarderQst=new float[QmobDielBoarderQstNum[2]];
		dielectricXQmobDBQst = new float[QmobDielBoarderQstNum[2]];
		dielectricYQmobDBQst = new float[QmobDielBoarderQstNum[2]];
		dielectricZQmobDBQst = new float[QmobDielBoarderQstNum[2]];
		dielectricXmQmobDBQst = new float[QmobDielBoarderQstNum[2]];
		dielectricYmQmobDBQst = new float[QmobDielBoarderQstNum[2]];
		dielectricZmQmobDBQst = new float[QmobDielBoarderQstNum[2]];
		dielectricZQmobDBSUMQst = new float[QmobDielBoarderQstNum[2]];
	}
	
	fprintf(stdout,"		QmobNum:............... [%d,%d,%d]\n", QmobNum[0], QmobNum[1], QmobNum[2]);
	fprintf(stdout,"		QmobDielBoarderNum:.... [%d,%d,%d]\n", QmobDielBoarderNum[0], QmobDielBoarderNum[1], QmobDielBoarderNum[2]);
	fprintf(stdout,"		QmobDielBoarderNumQst:. [%d,%d,%d]\n", QmobDielBoarderQstNum[0], QmobDielBoarderQstNum[1], QmobDielBoarderQstNum[2]);
	
	float *Eps=World->NIndexing->Eps;
	iQmob=0;
	iQmob2=QmobNum[1];
	iQmobDielBoarder=0;
	iQmobDielBoarder2=QmobDielBoarderNum[1];
	iQmobDielBoarderQst=0;
	iQmobDielBoarderQst2=QmobDielBoarderQstNum[1];
	int QCount=0;
	for(k=1;k<GS_Z-1;k++) 
	{
		kgrid = k*GS_XY;
		for(j=1;j<GS_Y-1;j++) 
		{
			jgrid = kgrid+j*GS_X;
			for(i=1;i<GS_X-1;i++) 
			{
				GrdPnt = jgrid+i;
				if(NIndex[GrdPnt]&specChargeMask)
					QCount++;
				
				if(CalcVolume!=NULL)
					bCalcVolume=CalcVolume[GrdPnt];
				
				bQmobHere=false;
				for(IType=0;IType<World->NIonsTypes;IType++)
					if(World->C[IType][GrdPnt]>0.0)bQmobHere=true;
				
				if(bCalcVolume&&bQmobHere)
				{
					if(NIndex[GrdPnt]&BlackAndWhiteMask)
					{
						if((NIndex[GrdPnt]&ChargeDielBoarderMask)==ChargeDielBoarderMask)
						{
							IndexQmobDielBoarderQst[iQmobDielBoarderQst]=GrdPnt;
							dielectricXQmobDBQst[iQmobDielBoarderQst]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
							dielectricYQmobDBQst[iQmobDielBoarderQst]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
							dielectricZQmobDBQst[iQmobDielBoarderQst]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
							dielectricXmQmobDBQst[iQmobDielBoarderQst]=Eps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
							dielectricYmQmobDBQst[iQmobDielBoarderQst]=Eps[(NIndex[GrdPnt-GS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
							dielectricZmQmobDBQst[iQmobDielBoarderQst]=Eps[(NIndex[GrdPnt-GS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
							dielectricZQmobDBSUMQst[iQmobDielBoarderQst]=dielectricXQmobDBQst[iQmobDielBoarderQst]+dielectricYQmobDBQst[iQmobDielBoarderQst]+dielectricZQmobDBQst[iQmobDielBoarderQst]+dielectricXmQmobDBQst[iQmobDielBoarderQst]+dielectricYmQmobDBQst[iQmobDielBoarderQst]+dielectricZmQmobDBQst[iQmobDielBoarderQst];
							dielectricZQmobDBSUMQst[iQmobDielBoarderQst]=Relaxation/dielectricZQmobDBSUMQst[iQmobDielBoarderQst];
							QstQmobDielBoarderQst[iQmobDielBoarderQst]=NIndexing->Q[QCount];
							iQmobDielBoarderQst++;
						}
						else if(NIndex[GrdPnt]&DielBoarderMask)
						{
							IndexQmobDielBoarder[iQmobDielBoarder]=GrdPnt;
							dielectricXQmobDB[iQmobDielBoarder]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
							dielectricYQmobDB[iQmobDielBoarder]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
							dielectricZQmobDB[iQmobDielBoarder]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
							dielectricXmQmobDB[iQmobDielBoarder]=Eps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
							dielectricYmQmobDB[iQmobDielBoarder]=Eps[(NIndex[GrdPnt-GS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
							dielectricZmQmobDB[iQmobDielBoarder]=Eps[(NIndex[GrdPnt-GS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
							dielectricZQmobDBSUM[iQmobDielBoarder]=dielectricXQmobDB[iQmobDielBoarder]+dielectricYQmobDB[iQmobDielBoarder]+dielectricZQmobDB[iQmobDielBoarder]+dielectricXmQmobDB[iQmobDielBoarder]+dielectricYmQmobDB[iQmobDielBoarder]+dielectricZmQmobDB[iQmobDielBoarder];
							dielectricZQmobDBSUM[iQmobDielBoarder]=Relaxation/dielectricZQmobDBSUM[iQmobDielBoarder];
							iQmobDielBoarder++;
						}
						else
						{
							IndexQmob[iQmob]=GrdPnt;
							dielectricChMob[iQmob]=World->NIndexing->GetDielFloat(0,GrdPnt);
							iQmob++;
						}
					}
					else
					{
						if((NIndex[GrdPnt]&ChargeDielBoarderMask)==ChargeDielBoarderMask)
						{
							IndexQmobDielBoarderQst[iQmobDielBoarderQst2]=GrdPnt;
							dielectricXQmobDBQst[iQmobDielBoarderQst2]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
							dielectricYQmobDBQst[iQmobDielBoarderQst2]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
							dielectricZQmobDBQst[iQmobDielBoarderQst2]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
							dielectricXmQmobDBQst[iQmobDielBoarderQst2]=Eps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
							dielectricYmQmobDBQst[iQmobDielBoarderQst2]=Eps[(NIndex[GrdPnt-GS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
							dielectricZmQmobDBQst[iQmobDielBoarderQst2]=Eps[(NIndex[GrdPnt-GS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
							dielectricZQmobDBSUMQst[iQmobDielBoarderQst2]=dielectricXQmobDBQst[iQmobDielBoarderQst2]+dielectricYQmobDBQst[iQmobDielBoarderQst2]+dielectricZQmobDBQst[iQmobDielBoarderQst2]+dielectricXmQmobDBQst[iQmobDielBoarderQst2]+dielectricYmQmobDBQst[iQmobDielBoarderQst2]+dielectricZmQmobDBQst[iQmobDielBoarderQst2];
							dielectricZQmobDBSUMQst[iQmobDielBoarderQst2]=Relaxation/dielectricZQmobDBSUMQst[iQmobDielBoarderQst2];
							QstQmobDielBoarderQst[iQmobDielBoarderQst2]=NIndexing->Q[QCount];
							iQmobDielBoarderQst2++;
						}
						else if(NIndex[GrdPnt]&DielBoarderMask)
						{
							IndexQmobDielBoarder[iQmobDielBoarder2]=GrdPnt;
							dielectricXQmobDB[iQmobDielBoarder2]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
							dielectricYQmobDB[iQmobDielBoarder2]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
							dielectricZQmobDB[iQmobDielBoarder2]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
							dielectricXmQmobDB[iQmobDielBoarder2]=Eps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
							dielectricYmQmobDB[iQmobDielBoarder2]=Eps[(NIndex[GrdPnt-GS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
							dielectricZmQmobDB[iQmobDielBoarder2]=Eps[(NIndex[GrdPnt-GS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
							dielectricZQmobDBSUM[iQmobDielBoarder2]=dielectricXQmobDB[iQmobDielBoarder2]+dielectricYQmobDB[iQmobDielBoarder2]+dielectricZQmobDB[iQmobDielBoarder2]+dielectricXmQmobDB[iQmobDielBoarder2]+dielectricYmQmobDB[iQmobDielBoarder2]+dielectricZmQmobDB[iQmobDielBoarder2];
							dielectricZQmobDBSUM[iQmobDielBoarder2]=Relaxation/dielectricZQmobDBSUM[iQmobDielBoarder2];
							iQmobDielBoarder2++;
						}
						else
						{
							IndexQmob[iQmob2]=GrdPnt;
							dielectricChMob[iQmob2]=World->NIndexing->GetDielFloat(0,GrdPnt);
							iQmob2++;
						}
					}
				}
			}
		}
	}
	SetQmobFromConcentration();
	return EXIT_SUCCESS;
}


int PoissonSolver::Solve()
{
	if(solver==Auto)
	{
		if(World->NIndexing==NULL)return PoissonSolverAD();
		else return PoissonSolverNIB();
	}
	else if(solver==NodeIndexBased)
		return PoissonSolverNIB();
	else if(solver==ArrayDirect)
		return PoissonSolverAD();
	else if(solver==PNPC)
		return PoissonSolverW();
	return EXIT_FAILURE;
}
int PoissonSolver::PoissonSolverAD()
{
	float gridScale;
	int gridSizeX;
	int gridSizeY;
	int gridSizeZ;
	int iteration;
	int i,j,k;
	int gridPoint;
	float om1,om2,om2d6;
	float * potential;
	float * staticCharge;
	int gridSizeXY;
	int gridSizeXYZ;
	float temp1,temp2,temp3,temp4,temp5,temp6,temp7;
	float dynamicChargeFactor,IonStrengthFactor;
	float fpoh;
	double totalEnergyOld=totalEnergy;
	bool *PeriodicBoundaryCondition;
	float *positiveCharge;
	float *negativeCharge;
	
	float * dielectric[3];
	
	
	for(i=0;i<3;i++)
		dielectric[i] = World->Epsilon[i];
	
	if(!(World->Potential)||!(World->Qstat)){
		fprintf(stderr,"ERROR 110: Arrays NOT yet initialize\n");
		exit(105);
	}

	//cout<<"TRACE: CPoisson::poissonSolver() START\n";
	gridScale = World->GridScale;
	gridSizeX = World->GridSize[0];
	gridSizeY = World->GridSize[1];
	gridSizeZ = World->GridSize[2];
	gridSizeXY = gridSizeX*gridSizeY;
	gridSizeXYZ = gridSizeXY*gridSizeZ;
	
	PeriodicBoundaryCondition=World->PBC;
	
	potential = World->Potential;
	if(World->C==NULL)
	{
		positiveCharge = NULL;
		negativeCharge = NULL;
	}
	else
	{
		positiveCharge = World->C[0];
		negativeCharge = World->C[1];
	}

	om2 = Relaxation;
	om1 = 1.0-om2;
	om2d6 = om2/6.0;
	
	//poissonBoltzmannSolver: assuming cation and anion concentration profiles are the same
	//poissonBoltzmannSolver: using cation concentration profile to calculate Debye lengths
	fpoh = 4.0*M_PI*gridScale;
	
	/*convert from charge density to M then convert from M to k'2
	* (or inverse debyle length squ					ared).
	*/
	dynamicChargeFactor = (float)1.0/(COANGS*4.0*M_PI*DFACT*DFACT);
	IonStrengthFactor = 1.0/(4.0*M_PI*gridScale);
	
	staticCharge=ChargeSum;
	
	//For speading calc denom first
	if(positiveCharge==0)
		for(k=0;k<gridSizeXYZ;k++){
		staticCharge[k]=World->Qstat[k];
		}
		else
			for(k=0;k<gridSizeXYZ;k++){
			staticCharge[k]=World->Qstat[k]+positiveCharge[k]-negativeCharge[k];
			}
			for(i=0;i<ChargeNum[2];i++){
				gridPoint=IndexCharge[i];
				staticCharge[gridPoint]/=dielectricCh[i];
			}
	
			for(iteration=1;iteration<=MaxIterations;iteration++) {
				for(j=0;j<=1;j++){
					for(i=SingularNum[j];i<SingularNum[j+1];i++){
						gridPoint=IndexSingular[i];
						temp1 = potential[gridPoint+1]*dielectricXS[i];
						temp2 = potential[gridPoint-1]*dielectricXmS[i];
						temp3 = potential[gridPoint+gridSizeX]*dielectricYS[i];
						temp4 = potential[gridPoint-gridSizeX]*dielectricYmS[i];
						temp5 = potential[gridPoint+gridSizeXY]*dielectricZS[i];
						temp6 = potential[gridPoint-gridSizeXY]*dielectricZmS[i];				
						temp1 = temp1+temp2;
						temp2 = temp3+temp4;
						temp3 = temp5+temp6;
						temp1 = temp1+temp2;
						temp2 = temp3+staticCharge[gridPoint];								
						temp1 = dielectricZSSUM[i]*(temp1+temp2);
						potential[gridPoint] = potential[gridPoint]*om1+temp1;
					}
					for(i=ChargeNum[j];i<ChargeNum[j+1];i++){
						gridPoint=IndexCharge[i];
						temp1 = potential[gridPoint+1]+potential[gridPoint-1];
						temp2 = potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX];
						temp3 = potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY];
						temp4 = potential[gridPoint]*om1;
						temp5 = om2d6*(temp1+temp2+temp3+staticCharge[gridPoint]);
				//temp6 = denominator[gridPoint]*(temp3+staticCharge[gridPoint]);
						potential[gridPoint] = temp4+temp5;
					}
					for(i=DielBoarderNum[j];i<DielBoarderNum[j+1];i++){
						gridPoint=IndexDielBoarder[i];
						temp1 = potential[gridPoint+1]*dielectricXDB[i];
						temp2 = potential[gridPoint-1]*dielectricXmDB[i];
						temp3 = potential[gridPoint+gridSizeX]*dielectricYDB[i];
						temp4 = potential[gridPoint-gridSizeX]*dielectricYmDB[i];
						temp5 = potential[gridPoint+gridSizeXY]*dielectricZDB[i];
						temp6 = potential[gridPoint-gridSizeXY]*dielectricZmDB[i];
						temp7 = potential[gridPoint]*om1;
						potential[gridPoint] = temp7+dielectricZDBSUM[i]*(temp1+temp2+temp3+temp4+temp5+temp6);
					}
					for(i=NoSingularNum[j];i<NoSingularNum[j+1];i++){
						gridPoint=IndexNoSingular[i];
						potential[gridPoint] = potential[gridPoint]*om1 + om2d6 * (potential[gridPoint+1] + potential[gridPoint-1] + potential[gridPoint+gridSizeX] + potential[gridPoint-gridSizeX] + potential[gridPoint+gridSizeXY] + potential[gridPoint-gridSizeXY]);
					}
					World->BorderExchange(potential);
				}
		
		//World->BorderExchange(potential);
		
				if((verbose&&(iteration%ConvergenceCheck==0))||iteration==MaxIterations)
				{
					temp1=totalEnergy;
					totalEnergy=CalculateEnergyPAD(fpoh,potential,staticCharge,dielectricCh, IndexCharge, IndexSingular, ChargeNum[2], SingularNum[2]);
					totalChange = fabs(totalEnergy-temp1);
					relativeChange=totalChange/totalEnergy;
					ConvFac=totalChange;
					
					if(verbose)
						pnpPrintGroup0("<PoissonIterations Nit=\"%8d\" E=\"%20.16le\" dE=\"%17.8lg\" rel.E=\"%17.8le\" ConvFac=\"%17.8le\"/>\n", iteration, totalEnergy, totalChange, relativeChange,ConvFac);

// #ifdef WIN32
//			 if(totalEnergy>1E13)
// #else
//			 if(totalEnergy>1E13||totalEnergy==NAN)
// #endif
					if(totalEnergy>1E13)
					{
						fprintf(stdout,"totalEnergy>1E13 || totalEnergy==NAN />\n", totalEnergy);
						return EXIT_FAILURE;
					}
					
					if(ConvFac<Convergence&&iteration>MinIterations)iteration = MaxIterations+1;
				}
			}
	//Calculate last energy
			totalEnergy=CalculateEnergyPAD(fpoh,potential,staticCharge,dielectricCh, IndexCharge, IndexSingular,ChargeNum[2], SingularNum[2]);
	
			if(verbose)
				pnpPrintGroup0("<PoissonFinal E=\"%.10e\" Eerr=\"%.10e\" Niter=\"%d\"/>\n", totalEnergy, totalChange, iteration-1);
			return EXIT_SUCCESS;
}
float PoissonSolver::CalculateEnergyPAD(float fpoh,float *Potential,float *StaticCharge,float *Epsilon,int *IndexCharge, int *IndexSingular,int ChargeNum,int SingularNum)
{
	int i,gridPoint;
	double Energy = 0.0;

	for(i=0;i<ChargeNum;i++){
		gridPoint=IndexCharge[i];
		Energy+= Potential[gridPoint]*StaticCharge[gridPoint]*Epsilon[i];
	}
	for(i=0;i<SingularNum;i++){
		gridPoint=IndexSingular[i];
		Energy+= Potential[gridPoint]*StaticCharge[gridPoint];
	}
#ifdef MPI_PARALLEL
	if(World->NProcs!=1)
	{
		MPI::Status	status;
		double EnergyProc;
		if(World->MyRank==0)
		{
			for(i=1;i<World->NProcs;i++)
			{
				pnpsapp->MyComGroup.Recv(&EnergyProc, 1, MPI::DOUBLE, i, SEND_ENERGY, status);
				Energy+=EnergyProc;
			}
		}
		else
		{
			pnpsapp->MyComGroup.Send(&Energy, 1, MPI::DOUBLE, 0, SEND_ENERGY);
		}
	}
#endif
	//Energy=Energy/(fpoh*2.0);
	World->SystemEnergy=Energy/(fpoh*2.0);
	return Energy/(fpoh*2.0);
}
int PoissonSolver::PoissonSolverW()
{
	int gridSizeX;
	int gridSizeY;
	int gridSizeZ;
	float gridScale;
	int iteration;
	int i,j,k;
	long gridPoint;
	float om1,om2,om2six;
	float * potential;
	float * potentialToo;
	float * dielectric[3];
	long gridSizeXY;
	long gridSizeXYZ;
	double totalChange;
	double totalEnergy;
	float change;
	long nextBorderPoint;
	long currentBorderCount;
	long * borderPoints;
	short * typeOfBorderPoint;
	float * om2InverseDielectricSum;
	float * positiveCharge;
	float * negativeCharge;
	float * chargeSum;
	long kgrid;
	long gridPointTemp;
	float temp1,temp2,temp3,temp4,temp5,temp6,temp7;
	int maxIterations;
	float convergence;
	long gridPointMax;
	float fpoh;

//	 assert(output!=NULL);
//	 assert(output->dielectricMap[0]!=NULL);
//	 assert(output->dielectricMap[1]!=NULL);
//	 assert(output->dielectricMap[2]!=NULL);
//	 assert(output->potentialMap!=NULL);

	assert(poissonSolverData!=NULL);
	assert(poissonSolverData->borderPoints!=NULL);

	gridScale = World->GridScale;
	gridSizeX = World->GridSize[0];
	gridSizeY = World->GridSize[1];
	gridSizeZ = World->GridSize[2];
	gridSizeXY = gridSizeX*gridSizeY;
	gridSizeXYZ = gridSizeXY*gridSizeZ;

	borderPoints = poissonSolverData->borderPoints;
	typeOfBorderPoint = poissonSolverData->typeOfBorderPoint;
	om2InverseDielectricSum = poissonSolverData->om2InverseDielectricSum;

	maxIterations = poissonSolverData->maxIterations;
	convergence = poissonSolverData->convergence;
	
	potential = World->Potential;
	potentialToo = potential;
	
	chargeSum = World->Qstat;
	if(World->C!=NULL)
	{
		positiveCharge = World->C[0];
		negativeCharge = World->C[1];
	}
	else
	{
		positiveCharge = NULL;
		negativeCharge = NULL;
		poissonSolverData->hasDynamicCharges=false;
	}
	
	for(i=0;i<3;i++)
		dielectric[i] = World->Epsilon[i];

	om2 = poissonSolverData->relaxation;
	om2six = om2/6;
	om1 = 1-om2;

	fpoh = 4*M_PI*gridScale;

	if(maxIterations>=POISSON_SOLVER_CONVERGENCE_CHECK) {
		fprintf(stdout,"poissonSolver: energy change (kT)		total energy (kT)		iteration\n");
		fprintf(stdout,"poissonSolver:-----------------------------------------------------\n");
	}

	if(poissonSolverData->hasDynamicCharges) {
		nextBorderPoint = borderPoints[0];
		currentBorderCount = 0;
		
		for(k=1;k<gridSizeZ-1;k++) {
			kgrid = k*gridSizeXY;
			for(j=1;j<gridSizeY-1;j++) {
				gridPoint = kgrid+j*gridSizeX;
				gridPointMax = gridPoint+gridSizeX-1;
				if(gridPointMax<nextBorderPoint) {
					for(gridPoint++;gridPoint<gridPointMax;gridPoint++) {
						chargeSum[gridPoint]+=positiveCharge[gridPoint]-negativeCharge[gridPoint];
						chargeSum[gridPoint]*=om2six/dielectric[0][gridPoint];
					}
				}
				else {
					for(gridPoint++;gridPoint<gridPointMax;gridPoint++) {
						if(gridPoint<nextBorderPoint) {
							chargeSum[gridPoint]+=positiveCharge[gridPoint]-negativeCharge[gridPoint];
							chargeSum[gridPoint]*=om2six/dielectric[0][gridPoint];
						}
						else {
							chargeSum[gridPoint]+=positiveCharge[gridPoint]-negativeCharge[gridPoint];
							currentBorderCount++;
							nextBorderPoint = borderPoints[currentBorderCount];
						}
					}
				}
			}
		}

		for(iteration=1;iteration<=MaxIterations;iteration++) {
			nextBorderPoint = borderPoints[0];
			currentBorderCount = 0;

			if(iteration%ConvergenceCheck==0) {
				totalChange = 0;
				totalEnergy = 0;
				for(k=1;k<gridSizeZ-1;k++) {
					kgrid = k*gridSizeXY;
					for(j=1;j<gridSizeY-1;j++) {
						gridPoint = kgrid+j*gridSizeX;
						gridPointMax = gridPoint+gridSizeX-1;
						while(nextBorderPoint<gridPointMax) {
							for(gridPoint++;gridPoint<nextBorderPoint;gridPoint++) {
								temp1 = potential[gridPoint+1]+potential[gridPoint-1];
								temp2 = potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX];
								temp3 = potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY];
								temp4 = potential[gridPoint]*om2;
								temp5 = om2six*(temp1+temp2);
								temp6 = om2six*temp3-temp4;
								change = temp5+temp6+chargeSum[gridPoint];
								potentialToo[gridPoint] = potential[gridPoint]+change;
								totalChange+=change*chargeSum[gridPoint]*dielectric[0][gridPoint]/om2six;
								totalEnergy+=potential[gridPoint]*chargeSum[gridPoint]*dielectric[0][gridPoint]/om2six;
							}
							temp1 = potential[gridPoint+1]*dielectric[0][gridPoint];
							temp2 = potential[gridPoint-1]*dielectric[0][gridPoint-1];
							temp3 = potential[gridPoint+gridSizeX]*dielectric[1][gridPoint];
							temp4 = potential[gridPoint-gridSizeX]*dielectric[1][gridPoint-gridSizeX];
							temp5 = potential[gridPoint+gridSizeXY]*dielectric[2][gridPoint];
							temp6 = potential[gridPoint-gridSizeXY]*dielectric[2][gridPoint-gridSizeXY];
							temp7 = potential[gridPoint]*om2;
							change = om2InverseDielectricSum[currentBorderCount]*(temp1+temp2+temp3+temp4+temp5+temp6+chargeSum[gridPoint])-temp7;
							potentialToo[gridPoint] = potential[gridPoint]+change;
							totalChange+=change*chargeSum[gridPoint];
							totalEnergy+=potential[gridPoint]*chargeSum[gridPoint];
							currentBorderCount++;
							nextBorderPoint = borderPoints[currentBorderCount];
						}
						for(gridPoint++;gridPoint<gridPointMax;gridPoint++) {
							temp1 = potential[gridPoint+1]+potential[gridPoint-1];
							temp2 = potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX];
							temp3 = potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY];
							temp4 = potential[gridPoint]*om2;
							temp5 = om2six*(temp1+temp2);
							temp6 = om2six*temp3-temp4;
							change = temp5+temp6+chargeSum[gridPoint];
							potentialToo[gridPoint] = potential[gridPoint]+change;
							totalChange+=change*chargeSum[gridPoint]*dielectric[0][gridPoint]/om2six;
							totalEnergy+=potential[gridPoint]*chargeSum[gridPoint]*dielectric[0][gridPoint]/om2six;
						}
					}
				}
	
				totalEnergy = totalEnergy/(fpoh*2);
				totalChange = fabs(totalChange/(fpoh*2));
				fprintf(stdout,"poissonSolver: %18e		%17e	 %10i\n",totalChange,totalEnergy,iteration);
	
				if(totalChange<Convergence)
					iteration = MaxIterations+1;
			}
			else {
				for(k=1;k<gridSizeZ-1;k++) {
					kgrid = k*gridSizeXY;
					for(j=1;j<gridSizeY-1;j++) {
						gridPoint = kgrid+j*gridSizeX;
						gridPointMax = gridPoint+gridSizeX-1;
						while(nextBorderPoint<gridPointMax) {
							gridPointTemp = gridPoint;
							for(gridPoint++;gridPoint<nextBorderPoint;gridPoint+=2) {
								potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY])+chargeSum[gridPoint];
							}
							for(gridPoint=gridPointTemp+2;gridPoint<nextBorderPoint;gridPoint+=2) {
								potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY])+chargeSum[gridPoint];
							}
							gridPoint = nextBorderPoint;
							temp1 = potential[gridPoint+1]*dielectric[0][gridPoint];
							temp2 = potential[gridPoint-1]*dielectric[0][gridPoint-1];
							temp3 = potential[gridPoint+gridSizeX]*dielectric[1][gridPoint];
							temp4 = potential[gridPoint-gridSizeX]*dielectric[1][gridPoint-gridSizeX];
							temp5 = potential[gridPoint+gridSizeXY]*dielectric[2][gridPoint];
							temp6 = potential[gridPoint-gridSizeXY]*dielectric[2][gridPoint-gridSizeXY];
							temp7 = potential[gridPoint]*om1;
							potentialToo[gridPoint] = temp7+om2InverseDielectricSum[currentBorderCount]*(temp1+temp2+temp3+temp4+temp5+temp6+chargeSum[gridPoint]);
							currentBorderCount++;
							nextBorderPoint = borderPoints[currentBorderCount];
						}
						gridPointTemp = gridPoint;
						for(gridPoint++;gridPoint<gridPointMax;gridPoint+=2) {
							potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY])+chargeSum[gridPoint];
						}
						for(gridPoint=gridPointTemp+2;gridPoint<gridPointMax;gridPoint+=2) {
							potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY])+chargeSum[gridPoint];
						}
					}
				}
			}
		}
		
		nextBorderPoint = borderPoints[0];
		currentBorderCount = 0;
		
		for(k=1;k<gridSizeZ-1;k++) {
			kgrid = k*gridSizeXY;
			for(j=1;j<gridSizeY-1;j++) {
				gridPoint = kgrid+j*gridSizeX;
				gridPointMax = gridPoint+gridSizeX-1;
				if(gridPointMax<nextBorderPoint) {
					for(gridPoint++;gridPoint<gridPointMax;gridPoint++) {
						chargeSum[gridPoint]/=om2six/dielectric[0][gridPoint];
						chargeSum[gridPoint]-=positiveCharge[gridPoint]-negativeCharge[gridPoint];
					}
				}
				else {
					for(gridPoint++;gridPoint<gridPointMax;gridPoint++) {
						if(gridPoint<nextBorderPoint) {
							chargeSum[gridPoint]/=om2six/dielectric[0][gridPoint];
							chargeSum[gridPoint]-=positiveCharge[gridPoint]-negativeCharge[gridPoint];
						}
						else {
							chargeSum[gridPoint]-=positiveCharge[gridPoint]-negativeCharge[gridPoint];
							currentBorderCount++;
							nextBorderPoint = borderPoints[currentBorderCount];
						}
					}
				}
			}
		}
	}
	else {
		currentBorderCount = 0;
		
		for(gridPoint=borderPoints[currentBorderCount];gridPoint<gridSizeXYZ;gridPoint=borderPoints[++currentBorderCount]) {
			if(typeOfBorderPoint[currentBorderCount]==POISSON_SOLVER_CHARGED_POINT) {
				chargeSum[gridPoint]/=dielectric[0][gridPoint];
			}
		}
			
		for(iteration=1;iteration<=maxIterations;iteration++) {
			nextBorderPoint = borderPoints[0];
			currentBorderCount = 0;

			if(iteration%ConvergenceCheck==0) {
				totalChange = 0;
				totalEnergy = 0;
				for(k=1;k<gridSizeZ-1;k++) {
					kgrid = k*gridSizeXY;
					for(j=1;j<gridSizeY-1;j++) {
						gridPoint = kgrid+j*gridSizeX;
						gridPointMax = gridPoint+gridSizeX-1;
						while(nextBorderPoint<gridPointMax) {
							gridPointTemp = gridPoint;
							for(gridPoint++;gridPoint<nextBorderPoint;gridPoint+=2) {
								potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY]);
							}
							for(gridPoint=gridPointTemp+2;gridPoint<nextBorderPoint;gridPoint+=2) {
								potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY]);
							}
							gridPoint = nextBorderPoint;
							if(typeOfBorderPoint[currentBorderCount]==POISSON_SOLVER_CHARGED_POINT) {
								temp1 = potential[gridPoint+1]+potential[gridPoint-1];
								temp2 = potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX];
								temp3 = potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY];
								temp4 = potential[gridPoint]*om2;
								temp5 = om2six*(temp1+temp2);
								temp6 = om2six*(temp3+chargeSum[gridPoint]);
								change = temp5+temp6-temp4;
								potentialToo[gridPoint] = potential[gridPoint]+change;
								totalChange+=change*chargeSum[gridPoint]*dielectric[0][gridPoint];
								totalEnergy+=potential[gridPoint]*chargeSum[gridPoint]*dielectric[0][gridPoint];
							}
							else if(typeOfBorderPoint[currentBorderCount]==POISSON_SOLVER_UNCOMMON_DIELECTRIC_POINT) {
								temp1 = potential[gridPoint+1]*dielectric[0][gridPoint];
								temp2 = potential[gridPoint-1]*dielectric[0][gridPoint-1];
								temp3 = potential[gridPoint+gridSizeX]*dielectric[1][gridPoint];
								temp4 = potential[gridPoint-gridSizeX]*dielectric[1][gridPoint-gridSizeX];
								temp5 = potential[gridPoint+gridSizeXY]*dielectric[2][gridPoint];
								temp6 = potential[gridPoint-gridSizeXY]*dielectric[2][gridPoint-gridSizeXY];
								temp7 = potential[gridPoint]*om1;
								potentialToo[gridPoint] = temp7+om2InverseDielectricSum[currentBorderCount]*(temp1+temp2+temp3+temp4+temp5+temp6);
							}
							else {
								temp1 = potential[gridPoint+1]*dielectric[0][gridPoint];
								temp2 = potential[gridPoint-1]*dielectric[0][gridPoint-1];
								temp3 = potential[gridPoint+gridSizeX]*dielectric[1][gridPoint];
								temp4 = potential[gridPoint-gridSizeX]*dielectric[1][gridPoint-gridSizeX];
								temp5 = potential[gridPoint+gridSizeXY]*dielectric[2][gridPoint];
								temp6 = potential[gridPoint-gridSizeXY]*dielectric[2][gridPoint-gridSizeXY];
								change = om2InverseDielectricSum[currentBorderCount]*(temp1+temp2+temp3+temp4+temp5+temp6+chargeSum[gridPoint])-potential[gridPoint]*om2;
								potentialToo[gridPoint] = potential[gridPoint]+change;
								totalChange+=change*chargeSum[gridPoint];
								totalEnergy+=potential[gridPoint]*chargeSum[gridPoint];
							}
							currentBorderCount++;
							nextBorderPoint = borderPoints[currentBorderCount];
						}
						gridPointTemp = gridPoint;
						for(gridPoint=gridPointTemp+1;gridPoint<gridPointMax;gridPoint+=2) {
							potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY]);
						}
						for(gridPoint=gridPointTemp+2;gridPoint<gridPointMax;gridPoint+=2) {
							potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY]);
						}
					}
				}
	
				totalEnergy = totalEnergy/(fpoh*2);
				totalChange = fabs(totalChange/(fpoh*2));
				relativeChange = totalChange/totalEnergy;
				
				if(verbose)
					pnpPrintGroup0("<PoissonIterations Nit=\"%8d\" E=\"%20.16le\" dE=\"%17.8lg\" rel.E=\"%17.8le\"/>\n", iteration, totalEnergy, totalChange, relativeChange);
				if(totalChange<convergence)
					iteration = maxIterations+1;
			}
			else {
				for(k=1;k<gridSizeZ-1;k++) {
					kgrid = k*gridSizeXY;
					for(j=1;j<gridSizeY-1;j++) {
						gridPoint = kgrid+j*gridSizeX;
						gridPointMax = gridPoint+gridSizeX-1;
						while(nextBorderPoint<gridPointMax) {
							gridPointTemp = gridPoint;
							for(gridPoint++;gridPoint<nextBorderPoint;gridPoint+=2) {
								potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY]);
							}
							for(gridPoint=gridPointTemp+2;gridPoint<nextBorderPoint;gridPoint+=2) {
								potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY]);
							}
							gridPoint = nextBorderPoint;
							if(typeOfBorderPoint[currentBorderCount]==POISSON_SOLVER_CHARGED_POINT) {
								temp1 = potential[gridPoint+1]+potential[gridPoint-1];
								temp2 = potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX];
								temp3 = potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY];
								temp4 = potential[gridPoint]*om1;
								temp5 = om2six*(temp1+temp2);
								temp6 = om2six*(temp3+chargeSum[gridPoint]);
								potentialToo[gridPoint] = temp4+temp5+temp6;
							}
							else if(typeOfBorderPoint[currentBorderCount]==POISSON_SOLVER_UNCOMMON_DIELECTRIC_POINT) {
								temp1 = potential[gridPoint+1]*dielectric[0][gridPoint];
								temp2 = potential[gridPoint-1]*dielectric[0][gridPoint-1];
								temp3 = potential[gridPoint+gridSizeX]*dielectric[1][gridPoint];
								temp4 = potential[gridPoint-gridSizeX]*dielectric[1][gridPoint-gridSizeX];
								temp5 = potential[gridPoint+gridSizeXY]*dielectric[2][gridPoint];
								temp6 = potential[gridPoint-gridSizeXY]*dielectric[2][gridPoint-gridSizeXY];
								temp7 = potential[gridPoint]*om1;
								potentialToo[gridPoint] = temp7+om2InverseDielectricSum[currentBorderCount]*(temp1+temp2+temp3+temp4+temp5+temp6);
							}
							else {
								temp1 = potential[gridPoint+1]*dielectric[0][gridPoint];
								temp2 = potential[gridPoint-1]*dielectric[0][gridPoint-1];
								temp3 = potential[gridPoint+gridSizeX]*dielectric[1][gridPoint];
								temp4 = potential[gridPoint-gridSizeX]*dielectric[1][gridPoint-gridSizeX];
								temp5 = potential[gridPoint+gridSizeXY]*dielectric[2][gridPoint];
								temp6 = potential[gridPoint-gridSizeXY]*dielectric[2][gridPoint-gridSizeXY];
								potentialToo[gridPoint] = potential[gridPoint]*om1+om2InverseDielectricSum[currentBorderCount]*(temp1+temp2+temp3+temp4+temp5+temp6+chargeSum[gridPoint]);
							}
							currentBorderCount++;
							nextBorderPoint = borderPoints[currentBorderCount];
						}
						gridPointTemp = gridPoint;
						for(gridPoint=gridPointTemp+1;gridPoint<gridPointMax;gridPoint+=2) {
							potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY]);
						}
						for(gridPoint=gridPointTemp+2;gridPoint<gridPointMax;gridPoint+=2) {
							potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY]);
						}
					}
				}
			}
		}

		currentBorderCount = 0;
		
		for(gridPoint=borderPoints[currentBorderCount];gridPoint<gridSizeXYZ;gridPoint=borderPoints[currentBorderCount++]) {
			if(typeOfBorderPoint[currentBorderCount]==POISSON_SOLVER_CHARGED_POINT) {
				chargeSum[gridPoint]*=dielectric[0][gridPoint];
			}
		}
	}


	return EXIT_SUCCESS;
}
int PoissonSolver::PoissonSolverNIB(bool ckenergy)
{
	//DbgPrint2("CPoisson::PoissonSolver %d\n",World->MyRank);
	//tmp var
	int iteration;
	int i,j,k;
	int GrdPnt;
	
#ifndef PNPDOUBLE
	float tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7;
#else
	double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7;
#endif
	

	//local vars
#ifndef PNPDOUBLE
	float gridScale = World->GridScale;
#else
	double gridScale = World->GridScale;
#endif
	int GS_X = World->GridSize[0];
	int GS_Y = World->GridSize[1];
	int GS_Z = World->GridSize[2];
	int GS_XY = GS_X*GS_Y;
	int GS_XYZ = GS_XY*GS_Z;
	bool *PeriodicBoundaryCondition = World->PBC;
#ifndef PNPDOUBLE
	//vars
	float fpoh = 4.0*M_PI*gridScale;
	float dynamicChargeFactor = 1.0/(COANGS*4.0*M_PI*DFACT*DFACT);
	float IonStrengthFactor = 1.0/(4.0*M_PI*gridScale);
	float om2 = Relaxation;
	float om1 = 1.0-om2;
	float om2d6 = om2/6.0;
	double totalEnergyOld=totalEnergy;
	//arrays
#else
	//vars
	double fpoh = 4.0*M_PI*gridScale;
	double dynamicChargeFactor = 1.0/(COANGS*4.0*M_PI*DFACT*DFACT);
	double IonStrengthFactor = 1.0/(4.0*M_PI*gridScale);
	double om2 = Relaxation;
	double om1 = 1.0-om2;
	double om2d6 = om2/6.0;
	
	double totalEnergyOld=totalEnergy;
	//arrays
	if(potential == NULL)
	{
		potential = new double[GS_XYZ];
		for(i=0;i<GS_XYZ;i++)
			potential[i]= (double)World->Potential[i];
	}
#endif

	
	//Check out is everything 
	if(!(World->Potential)){
		fprintf(stderr,"ERROR 110: Arrays NOT yet initialize\n");
		exit(105);
	}
	int countConvFacHistory=0;
	double *ConvFacHistory=NULL;
	if(ConvFacMaxHistory>1)
	{
		ConvFacHistory=new double[ConvFacMaxHistory];
		for(i=0;i<ConvFacMaxHistory;i++)
			ConvFacHistory[i]=1e10;
	}
	//Iteration itself
	for(iteration=1;iteration<=MaxIterations;iteration++)
	{
		//calculation over black and white nodes
		for(j=0;j<=1;j++){
			for(i=SingularNum[j];i<SingularNum[j+1];i++){
				GrdPnt=IndexSingular[i];
				tmp1 = potential[GrdPnt+1]*dielectricXS[i];
				tmp2 = potential[GrdPnt-1]*dielectricXmS[i];
				tmp3 = potential[GrdPnt+GS_X]*dielectricYS[i];
				tmp4 = potential[GrdPnt-GS_X]*dielectricYmS[i];
				tmp5 = potential[GrdPnt+GS_XY]*dielectricZS[i];
				tmp6 = potential[GrdPnt-GS_XY]*dielectricZmS[i];				
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3+tmp4;
				tmp3 = tmp5+tmp6;
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3+QstS[i];								
				tmp1 = dielectricZSSUM[i]*(tmp1+tmp2);
				potential[GrdPnt] = potential[GrdPnt]*om1+tmp1;
			}
			for(i=ChargeNum[j];i<ChargeNum[j+1];i++){
				GrdPnt=IndexCharge[i];
				tmp1 = potential[GrdPnt+1]+potential[GrdPnt-1];
				tmp2 = potential[GrdPnt+GS_X]+potential[GrdPnt-GS_X];
				tmp3 = potential[GrdPnt+GS_XY]+potential[GrdPnt-GS_XY];
				tmp4 = potential[GrdPnt]*om1;
				tmp5 = om2d6*(tmp1+tmp2+tmp3+Qst[i]);
				//tmp6 = denominator[GrdPnt]*(tmp3+staticCharge[GrdPnt]);
				potential[GrdPnt] = tmp4+tmp5;
			}
			for(i=DielBoarderNum[j];i<DielBoarderNum[j+1];i++){
				GrdPnt=IndexDielBoarder[i];
				tmp1 = potential[GrdPnt+1]*dielectricXDB[i];
				tmp2 = potential[GrdPnt-1]*dielectricXmDB[i];
				tmp3 = potential[GrdPnt+GS_X]*dielectricYDB[i];
				tmp4 = potential[GrdPnt-GS_X]*dielectricYmDB[i];
				tmp5 = potential[GrdPnt+GS_XY]*dielectricZDB[i];
				tmp6 = potential[GrdPnt-GS_XY]*dielectricZmDB[i];
				tmp7 = potential[GrdPnt]*om1;
				potential[GrdPnt] = tmp7+dielectricZDBSUM[i]*(tmp1+tmp2+tmp3+tmp4+tmp5+tmp6);
				//potential[GrdPnt]+=100.0;
			}
			for(i=NoSingularNum[j];i<NoSingularNum[j+1];i++){
				GrdPnt=IndexNoSingular[i];
				potential[GrdPnt] = potential[GrdPnt]*om1 + om2d6 * (potential[GrdPnt+1] + potential[GrdPnt-1] + potential[GrdPnt+GS_X] + potential[GrdPnt-GS_X] + potential[GrdPnt+GS_XY] + potential[GrdPnt-GS_XY]);
				//potential[GrdPnt]+=10.0;
			}
			for(i=QmobNum[j];i<QmobNum[j+1];i++){
				GrdPnt=IndexQmob[i];
				tmp1 = potential[GrdPnt+1]+potential[GrdPnt-1];
				tmp2 = potential[GrdPnt+GS_X]+potential[GrdPnt-GS_X];
				tmp3 = potential[GrdPnt+GS_XY]+potential[GrdPnt-GS_XY];
				tmp4 = potential[GrdPnt]*om1;
				tmp5 = om2d6*(tmp1+tmp2+tmp3+Qmob[i]);
				//tmp6 = denominator[GrdPnt]*(tmp3+staticCharge[GrdPnt]);
				potential[GrdPnt] = tmp4+tmp5;
				//potential[GrdPnt]+=1.0;
			}
			for(i=QmobDielBoarderNum[j];i<QmobDielBoarderNum[j+1];i++){
				GrdPnt=IndexQmobDielBoarder[i];
				tmp1 = potential[GrdPnt+1]*dielectricXQmobDB[i];
				tmp2 = potential[GrdPnt-1]*dielectricXmQmobDB[i];
				tmp3 = potential[GrdPnt+GS_X]*dielectricYQmobDB[i];
				tmp4 = potential[GrdPnt-GS_X]*dielectricYmQmobDB[i];
				tmp5 = potential[GrdPnt+GS_XY]*dielectricZQmobDB[i];
				tmp6 = potential[GrdPnt-GS_XY]*dielectricZmQmobDB[i];				
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3+tmp4;
				tmp3 = tmp5+tmp6;
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3+QmobDielBoarder[i];								
				tmp1 = dielectricZQmobDBSUM[i]*(tmp1+tmp2);
				potential[GrdPnt] = potential[GrdPnt]*om1+tmp1;
			}
			for(i=QmobDielBoarderQstNum[j];i<QmobDielBoarderQstNum[j+1];i++){
				GrdPnt=IndexQmobDielBoarderQst[i];
				tmp1 = potential[GrdPnt+1]*dielectricXQmobDBQst[i];
				tmp2 = potential[GrdPnt-1]*dielectricXmQmobDBQst[i];
				tmp3 = potential[GrdPnt+GS_X]*dielectricYQmobDBQst[i];
				tmp4 = potential[GrdPnt-GS_X]*dielectricYmQmobDBQst[i];
				tmp5 = potential[GrdPnt+GS_XY]*dielectricZQmobDBQst[i];
				tmp6 = potential[GrdPnt-GS_XY]*dielectricZmQmobDBQst[i];				
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3+tmp4;
				tmp3 = tmp5+tmp6;
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3+QmobDielBoarderQst[i]+QstQmobDielBoarderQst[i];								
				tmp1 = dielectricZQmobDBSUMQst[i]*(tmp1+tmp2);
				potential[GrdPnt] = potential[GrdPnt]*om1+tmp1;
			}
#ifndef PNPDOUBLE
			World->BorderExchange(potential);
#else
			World->BorderExchangeDouble(potential);
#endif
		}
		//checking and printing energy
		
		if(ckenergy)
			if((verbose&&(iteration%ConvergenceCheck==0))||iteration==MaxIterations)
		{
			CalcSystemEnergy(iteration);
			relativeChange=totalChange/totalEnergy;
			
			if(verbose)
			{
				//pnpPrintGroup0("<PoissonIterations Nit=\"%8d\" E=\"%20.16e\" dE=\"%.4e\" rel.E=\"%.4e\" ConvFac=\"%.4e\"/>\n", iteration, totalEnergy, totalChange, relativeChange,ConvFac);
				if(iteration/ConvergenceCheck<=1)
				{
					pnpPrintGroup0("P     =========================================================================\n");
					pnpPrintGroup0("P      %9s %22s %12s %12s %12s\n","Iteration", "Energy,kT","dE","rel.E","ConvFac");
					pnpPrintGroup0("P     -------------------------------------------------------------------------\n");
				}
				pnpPrintGroup0("P      %9d %22.14e %12.4e %12.4e %12.4e\n", iteration, totalEnergy, totalChange, relativeChange,ConvFac);
			}
			if(totalEnergy>1E13)
			{
				pnpPrintGroup0("P     -------------------------------------------------------------------------\n");
				pnpPrintGroup0("P      ERROR: Poisson Solver has diverged, try smaller relaxation\n");
				pnpPrintGroup0("P     =========================================================================\n");
				return EXIT_FAILURE;
			}
			if(ConvFacMaxHistory>1)
			{
				ConvFacHistory[countConvFacHistory]=ConvFac;
				countConvFacHistory++;
				if(countConvFacHistory>=ConvFacMaxHistory)
					countConvFacHistory=0;
				
				bool convereged=true;
				for(i=0;i<ConvFacMaxHistory;i++)
					convereged=convereged&&(ConvFacHistory[i]<=Convergence);
				if(convereged)break;//iteration = MaxIterations+1;
			}
			else if(ConvFac<Convergence&&iteration>MinIterations)break;
		}
	}
	if(verbose&&ckenergy)
	{
		pnpPrintGroup0("P     -------------------------------------------------------------------------\n");
		pnpPrintGroup0("P      Results: E=%.14e Niter=%d\n", totalEnergy, iteration-1);
		pnpPrintGroup0("P     =========================================================================\n");
		//pnpPrintGroup0("<PoissonFinal E=\"%.10e\" Eerr=\"%.10e\" Niter=\"%d\"/>\n", totalEnergy, totalChange, iteration-1);
	}
#ifdef PNPDOUBLE
	for(i=0;i<GS_XYZ;i++)
		World->Potential[i] = (float)potential[i];
	delete [] potential;
#endif
	
	DeleteCArray(ConvFacHistory);
	
	if(Convergence>0.0&&ConvFac>Convergence)
		return EXIT_FAILURE;
	return EXIT_SUCCESS;
}
int PoissonSolver::CalcSystemEnergy(int iteration)
{
// 	CalcSystemEnergyFloat(iteration);
	if(WayToCalcSystemEnergy==0)
	{
		if(PMFCalculation0!=NULL)
			CalcSystemEnergy4PMF(iteration);
		else
			CalcSystemEnergyStdDevPhi(iteration);
	}
	else
	{
		if(WayToCalcSystemEnergy==1)
			CalcSystemEnergyStdDevPhi(iteration);
	}
	//CalcSystemEnergyDouble(iteration);
		//CalcSystemEnergyMaxPhiChange(iteration);
// 	CalcSystemEnergyLongDouble(iteration);
// 	CalcSystemEnergy0(iteration);
// 	CalcSystemEnergy1(iteration);
 	//CalcSystemEnergy3(iteration);
	//if(iteration==MaxIterations)
	//	CalcSystemEnergyAnalizer(iteration);
	return EXIT_SUCCESS;
}
int PoissonSolver::CalcSystemEnergyDouble(int iteration)
{
	int i,GrdPnt;
	static double oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	double Dev=0.0;
	double tmp,SumSQ=0.0,SumAbs=0.0;
	double EnergyCharge = 0.0;
	double EnergySingular = 0.0;
	double EnergyQmob = 0.0;
	double OldTotalEnergy=totalEnergy;
	float fpoh = 4.0*M_PI*World->GridScale;
	
	for(i=0;i<ChargeNum[2];i++)
	{
		GrdPnt=IndexCharge[i];
		
		tmp=double(potential[GrdPnt])*double(Qst[i])*double(dielectricCh[i]);
		EnergyCharge+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		tmp=double(QstS[i])*double(potential[GrdPnt]);
		EnergySingular+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobNum[2];i++)
	{
		GrdPnt=IndexQmob[i];
		tmp=double(potential[GrdPnt])*double(Qmob[i])*double(dielectricChMob[i]);
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobDielBoarderNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarder[i];
		tmp=double(potential[GrdPnt])*double(QmobDielBoarder[i]);
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	DbgPrint0("CSE_Double  :E=%.16e Eq=%.10e Esing=%.10e EQmob=%.10e\n"
			,(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));
#ifdef MPI_PARALLEL
	int dest;
	double locEnergyCharge;
	double locEnergySingular;
	double locEnergyQmob;
	double locSumAbs,locSumSQ;
	if(World->MyRank==0)
	{
		for(dest=1;dest<World->NProcs;dest++)
		{
			pnpsapp->MyComGroup.Recv(&locEnergyCharge, 1, MPI::DOUBLE, dest, 1);
			pnpsapp->MyComGroup.Recv(&locEnergySingular, 1, MPI::DOUBLE, dest, 2);
			pnpsapp->MyComGroup.Recv(&locEnergyQmob, 1, MPI::DOUBLE, dest, 3);
			
			
			EnergyCharge+=locEnergyCharge;
			EnergySingular+=locEnergySingular;
			EnergyQmob+=locEnergyQmob;
			pnpsapp->MyComGroup.Recv(&locSumAbs, 1, MPI::DOUBLE, dest, 4);
			pnpsapp->MyComGroup.Recv(&locSumSQ, 1, MPI::DOUBLE, dest, 5);
			SumSQ+=locSumSQ;
			SumAbs+=locSumAbs;
		}
		DbgPrint0("CSE_DoubleF :E=%.16e Eq=%.10e Esing=%.10e EQmob=%.10e\n"
			,(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));
	}
	else
	{
		pnpsapp->MyComGroup.Send(&EnergyCharge, 1, MPI::DOUBLE, 0, 1);
		pnpsapp->MyComGroup.Send(&EnergySingular, 1, MPI::DOUBLE, 0, 2);
		pnpsapp->MyComGroup.Send(&EnergyQmob, 1, MPI::DOUBLE, 0, 3);
		pnpsapp->MyComGroup.Send(&SumSQ, 1, MPI::DOUBLE, 0, 4);
		pnpsapp->MyComGroup.Send(&SumAbs, 1, MPI::DOUBLE, 0, 5);
	}
	
	pnpsapp->MyComGroup.Bcast(&EnergyCharge, 1, MPI::DOUBLE, 0);
	pnpsapp->MyComGroup.Bcast(&EnergySingular, 1, MPI::DOUBLE, 0);
	pnpsapp->MyComGroup.Bcast(&EnergyQmob, 1, MPI::DOUBLE, 0);
	pnpsapp->MyComGroup.Bcast(&SumSQ, 1, MPI::DOUBLE, 0);
	pnpsapp->MyComGroup.Bcast(&SumAbs, 1, MPI::DOUBLE, 0);
	
#endif

	totalEnergy=(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0);
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac = sqrt(fabs(oldSumSQ-SumSQ))/(fpoh*2.0);
	Dev=oldSumAbs-SumAbs;
	//if(oldSumSQ==0.0)fprintf(stdout,"PT Nit E dE RMSD dEabs\n");
	//if(verbose)fprintf(stdout,"PT %8d %17.12lg %17.12lg %17.12lg %17.12lg\n", iteration, totalEnergy, totalEnergy-oldTotEn,totalChange, Dev);
	
	oldSumSQ=SumSQ;
	oldSumAbs=SumAbs;
	oldTotEn=totalEnergy;
	World->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
int PoissonSolver::CalcSystemEnergyLongDouble(int iteration)
{
	int i,GrdPnt;
	static long double oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	long double Dev=0.0;
	long double tmp,SumSQ=0.0,SumAbs=0.0;
	long double EnergyCharge = 0.0;
	long double EnergySingular = 0.0;
	long double EnergyQmob = 0.0;
	long double OldTotalEnergy=totalEnergy;
	long double  fpoh = 4.0*M_PI*World->GridScale;
	
	for(i=0;i<ChargeNum[2];i++)
	{
		GrdPnt=IndexCharge[i];
		tmp=potential[GrdPnt];
		tmp=Qst[i];
		tmp=dielectricCh[i];
		
		tmp=potential[GrdPnt]*Qst[i]*dielectricCh[i];
		EnergyCharge+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		tmp=QstS[i]*potential[GrdPnt];
		EnergySingular+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobNum[2];i++)
	{
		GrdPnt=IndexQmob[i];
		tmp=potential[GrdPnt]*Qmob[i]*dielectricChMob[i];
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobDielBoarderNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarder[i];
		tmp=potential[GrdPnt]*QmobDielBoarder[i];
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	DbgPrint0("CSE_LDouble :E=%.16Le Eq=%.10Le Esing=%.10Le EQmob=%.10lle\n",(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));

	totalEnergy=(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0);
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac = sqrt(fabs(oldSumSQ-SumSQ))/(fpoh*2.0);
	Dev=oldSumAbs-SumAbs;
	//if(oldSumSQ==0.0)fprintf(stdout,"PT Nit E dE RMSD dEabs\n");
	//if(verbose)fprintf(stdout,"PT %8d %17.12lg %17.12lg %17.12lg %17.12lg\n", iteration, totalEnergy, totalEnergy-oldTotEn,totalChange, Dev);
	
	oldSumSQ=SumSQ;
	oldSumAbs=SumAbs;
	oldTotEn=totalEnergy;
	World->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
int PoissonSolver::CalcSystemEnergyFloat(int iteration)
{
	int i,GrdPnt;
	static float oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	float Dev=0.0;
	float tmp,SumSQ=0.0,SumAbs=0.0;
	float EnergyCharge = 0.0;
	float EnergySingular = 0.0;
	float EnergyQmob = 0.0;
	float OldTotalEnergy=totalEnergy;
	float  fpoh = 4.0*M_PI*World->GridScale;
	
	for(i=0;i<ChargeNum[2];i++)
	{
		GrdPnt=IndexCharge[i];
		tmp=potential[GrdPnt];
		tmp=Qst[i];
		tmp=dielectricCh[i];
		
		tmp=potential[GrdPnt]*Qst[i]*dielectricCh[i];
		EnergyCharge+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		tmp=QstS[i]*potential[GrdPnt];
		EnergySingular+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobNum[2];i++)
	{
		GrdPnt=IndexQmob[i];
		tmp=potential[GrdPnt]*Qmob[i]*dielectricChMob[i];
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobDielBoarderNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarder[i];
		tmp=potential[GrdPnt]*QmobDielBoarder[i];
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	
	DbgPrint0("CSE_Float   :E=%.16e Eq=%.10e Esing=%.10e EQmob=%.10e\n"
			,float((EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0)),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));

	totalEnergy=(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0);
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac = sqrt(fabs(oldSumSQ-SumSQ))/(fpoh*2.0);
	Dev=oldSumAbs-SumAbs;
	//if(oldSumSQ==0.0)fprintf(stdout,"PT Nit E dE RMSD dEabs\n");
	//if(verbose)fprintf(stdout,"PT %8d %17.12lg %17.12lg %17.12lg %17.12lg\n", iteration, totalEnergy, totalEnergy-oldTotEn,totalChange, Dev);
	
	oldSumSQ=SumSQ;
	oldSumAbs=SumAbs;
	oldTotEn=totalEnergy;
	World->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
int PoissonSolver::CalcSystemEnergy4PMF(int iteration)
{
	//CalcSystemEnergyDouble(iteration);
	int i,GrdPnt;
	static double oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	double Dev=0.0;
	double tmp,SumSQ=0.0,SumAbs=0.0;
	double EnergyCharge = 0.0;
	double EnergyMinusGsub = 0.0;
	double EnergyQmob = 0.0;
	double OldTotalEnergy=totalEnergy;
	double dphi;
	double  fpoh = 4.0*M_PI*World->GridScale;
	
	PNP_EXIT_FAIL_NULL(PMFCalculation0,"PMFCalculation not initialize\n");
	
	
	if(PMFCalculation0->Gsubtract!=NULL)
	{
		for(i=0;i<World->NIndexing->QNum;i++)
		{
			GrdPnt=PMFCalculation0->grdpntQst[i];
			tmp=World->NIndexing->Q[i]*potential[GrdPnt]+PMFCalculation0->Gsubtract[i];
			EnergyCharge+=tmp;
			//SumSQ+=tmp*tmp;
		}
	}
	else
	{
		for(i=0;i<World->NIndexing->QNum;i++)
		{
			GrdPnt=PMFCalculation0->grdpntQst[i];
			tmp=World->NIndexing->Q[i]*potential[GrdPnt];
			EnergyCharge+=tmp;
			//SumSQ+=tmp*tmp;
		}
	}
	SumSQ=0.0;
	for(i=0;i<ChargeNum[2];i++)
	{
		GrdPnt=IndexCharge[i];
		tmp=double(potential[GrdPnt])*double(Qst[i])*double(dielectricCh[i]);
		dphi=potential[GrdPnt]-PhiCharge[i];
		SumSQ+=dphi*dphi;
		PhiCharge[i]=potential[GrdPnt];
	}
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		tmp=double(QstS[i])*double(potential[GrdPnt]);
		dphi=potential[GrdPnt]-PhiSingular[i];
		SumSQ+=dphi*dphi;
		PhiSingular[i]=potential[GrdPnt];
	}
	//DbgPrint0("CSE_4PMF   :Eq=%.16e Em=%.16e\n"
	//		,float((EnergyCharge)/(fpoh*2.0)),(EnergyMinusGsub)/(fpoh*2.0));

	totalEnergy=(EnergyCharge)/(fpoh*2.0);
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac=sqrt(SumSQ/double(ChargeNum[2]+SingularNum[2]));
	//ConvFac = sqrt(fabs(oldSumSQ-SumSQ))/(fpoh*2.0);
	//Dev=oldSumAbs-SumAbs;
	//if(oldSumSQ==0.0)fprintf(stdout,"PT Nit E dE RMSD dEabs\n");
	//if(verbose)fprintf(stdout,"PT %8d %17.12lg %17.12lg %17.12lg %17.12lg\n", iteration, totalEnergy, totalEnergy-oldTotEn,totalChange, Dev);
	
	oldSumSQ=SumSQ;
	oldSumAbs=SumAbs;
	oldTotEn=totalEnergy;
	World->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
int PoissonSolver::CalcSystemEnergyMaxPhiChange(int iteration)
{
	int i,GrdPnt;
	static double oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	double Dev=0.0;
	double tmp,SumSQ=0.0,SumAbs=0.0;
	double EnergyCharge = 0.0;
	double EnergySingular = 0.0;
	double EnergyQmob = 0.0;
	double OldTotalEnergy=totalEnergy;
	double  fpoh = 4.0*M_PI*World->GridScale;
	float dphi=0.0,maxdphi=0.0;
	int GrdPntMaxDPhi=-1;
	
	for(i=0;i<ChargeNum[2];i++)
	{
		GrdPnt=IndexCharge[i];
		
		tmp=double(potential[GrdPnt])*double(Qst[i])*double(dielectricCh[i]);
		EnergyCharge+=tmp;
		dphi=fabs(potential[GrdPnt]-PhiCharge[i]);
		PhiCharge[i]=potential[GrdPnt];
		if(dphi>maxdphi)
		{
			maxdphi=dphi;
			GrdPntMaxDPhi=GrdPnt;
		}
	}
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		tmp=double(QstS[i])*double(potential[GrdPnt]);
		EnergySingular+=tmp;
		dphi=fabs(potential[GrdPnt]-PhiSingular[i]);
		PhiSingular[i]=potential[GrdPnt];
		if(dphi>maxdphi)
		{
			maxdphi=dphi;
			GrdPntMaxDPhi=GrdPnt;
		}
	}

	//DbgPrint0("CSE_0       :E=%.16e Eq=%.10e Esing=%.10e EQmob=%.10e\n",(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));
	DbgPrint0("CSE_MaxPhiChange %.5e at %d\n",maxdphi,GrdPntMaxDPhi);
	totalEnergy=(EnergyCharge+EnergySingular)/(fpoh*2.0);
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac=maxdphi;
	World->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
int PoissonSolver::CalcSystemEnergyStdDevPhi(int iteration)
{
	int i,GrdPnt;
	static double oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	double Dev=0.0;
	double tmp,SumSQ=0.0,SumAbs=0.0;
	double EnergyCharge = 0.0;
	double EnergySingular = 0.0;
	double EnergyQmob = 0.0;
	double OldTotalEnergy=totalEnergy;
	double  fpoh = 4.0*M_PI*World->GridScale;
	float dphi=0.0,maxdphi=0.0;
	int GrdPntMaxDPhi=-1;
	
	for(i=0;i<ChargeNum[2];i++)
	{
		GrdPnt=IndexCharge[i];
		
		tmp=double(potential[GrdPnt])*double(Qst[i])*double(dielectricCh[i]);
		EnergyCharge+=tmp;
		dphi=potential[GrdPnt]-PhiCharge[i];
		SumSQ+=dphi*dphi;
		PhiCharge[i]=potential[GrdPnt];
	}
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		tmp=double(QstS[i])*double(potential[GrdPnt]);
		EnergySingular+=tmp;
		dphi=potential[GrdPnt]-PhiSingular[i];
		SumSQ+=dphi*dphi;
		PhiSingular[i]=potential[GrdPnt];
	}

	//DbgPrint0("CSE_0       :E=%.16e Eq=%.10e Esing=%.10e EQmob=%.10e\n",(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));
	//DbgPrint0("CSE_MaxPhiChange %.5e at %d\n",maxdphi,GrdPntMaxDPhi);
	totalEnergy=(EnergyCharge+EnergySingular)/(fpoh*2.0);
	int ChargedNodes=ChargeNum[2]+SingularNum[2];
#ifdef MPI_PARALLEL
	int dest;
	pnpsapp->MyComGroup.Barrier();
			
	if(World->MyRank==0)
	{
		for(dest=1;dest<World->NProcs;dest++)
		{
			pnpsapp->MyComGroup.Recv(&tmp, 1, MPI::DOUBLE, dest, 0);
			totalEnergy+=tmp;
			pnpsapp->MyComGroup.Recv(&tmp, 1, MPI::DOUBLE, dest, 0);
			totalEnergy+=tmp;
			pnpsapp->MyComGroup.Recv(&i, 1, MPI::INT, dest, 0);
			ChargedNodes+=i;
		}
	}
	else
	{
		pnpsapp->MyComGroup.Send(&totalEnergy, 1, MPI::DOUBLE, 0, 0);
		pnpsapp->MyComGroup.Send(&SumSQ, 1, MPI::DOUBLE, 0, 0);
		pnpsapp->MyComGroup.Send(&ChargedNodes, 1, MPI::INT, 0, 0);
	}
#endif
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac=sqrt(SumSQ/double(ChargedNodes));
#ifdef MPI_PARALLEL
	//int dest;
	pnpsapp->MyComGroup.Barrier();
			
	if(World->MyRank==0)
	{
		for(dest=1;dest<World->NProcs;dest++)
		{
			pnpsapp->MyComGroup.Send(&totalEnergy, 1, MPI::DOUBLE, dest, 0);
			pnpsapp->MyComGroup.Send(&SumSQ, 1, MPI::DOUBLE, dest, 0);
			pnpsapp->MyComGroup.Send(&ChargedNodes, 1, MPI::INT, dest, 0);
			pnpsapp->MyComGroup.Send(&ConvFac, 1, MPI::DOUBLE, dest, 0);
		}
	}
	else
	{
		pnpsapp->MyComGroup.Recv(&totalEnergy, 1, MPI::DOUBLE, 0, 0);
		pnpsapp->MyComGroup.Recv(&SumSQ, 1, MPI::DOUBLE, 0, 0);
		pnpsapp->MyComGroup.Recv(&ChargedNodes, 1, MPI::INT, 0, 0);
		pnpsapp->MyComGroup.Recv(&ConvFac, 1, MPI::DOUBLE, 0, 0);
	}
#endif
	World->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
int PoissonSolver::CalcSystemEnergy0(int iteration)
{
	int i,GrdPnt;
	static double oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	double Dev=0.0;
	double tmp,SumSQ=0.0,SumAbs=0.0;
	double EnergyCharge = 0.0;
	double EnergySingular = 0.0;
	double EnergyQmob = 0.0;
	double OldTotalEnergy=totalEnergy;
	double  fpoh = 4.0*M_PI*World->GridScale;
	
	for(i=0;i<ChargeNum[2];i++)
	{
		GrdPnt=IndexCharge[i];
		
		tmp=double(potential[GrdPnt])*double(Qst[i])*double(dielectricCh[i]);
		EnergyCharge+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		tmp=double(QstS[i])*double(potential[GrdPnt]);
		EnergySingular+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobNum[2];i++)
	{
		GrdPnt=IndexQmob[i];
		tmp=double(potential[GrdPnt])*double(Qmob[i])*double(dielectricChMob[i]);
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobDielBoarderNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarder[i];
		tmp=double(potential[GrdPnt])*double(QmobDielBoarder[i]);
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	DbgPrint0("CSE_0       :E=%.16e Eq=%.10e Esing=%.10e EQmob=%.10e\n",(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));

	totalEnergy=(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0);
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac = sqrt(fabs(oldSumSQ-SumSQ))/(fpoh*2.0);
	Dev=oldSumAbs-SumAbs;
	//if(oldSumSQ==0.0)fprintf(stdout,"PT Nit E dE RMSD dEabs\n");
	//if(verbose)fprintf(stdout,"PT %8d %17.12lg %17.12lg %17.12lg %17.12lg\n", iteration, totalEnergy, totalEnergy-oldTotEn,totalChange, Dev);
	
	oldSumSQ=SumSQ;
	oldSumAbs=SumAbs;
	oldTotEn=totalEnergy;
	World->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
int PoissonSolver::CalcSystemEnergy3(int iteration)
{
	int i,GrdPnt;
	int AvrOver=10;
	static int Counter=0;
	static double oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	double Dev=0.0;
	double tmp,tmp2,SumSQ=0.0,SumAbs=0.0;
	double EnergyCharge = 0.0;
	double EnergySingular = 0.0;
	double EnergyQmob = 0.0;
	double OldTotalEnergy=totalEnergy;
	double  fpoh = 4.0*M_PI*World->GridScale;
	
	
	if(ChargesEnergy==NULL)
	{
		ChargesEnergy=new double[AvrOver];
		for(i=0;i<AvrOver;i++)ChargesEnergy[i]=-1e300;
		Counter=0;
	}
	
	for(i=0;i<ChargeNum[2];i++)
	{
		GrdPnt=IndexCharge[i];
		
		tmp=double(potential[GrdPnt])*double(Qst[i])*double(dielectricCh[i]);
		EnergyCharge+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		tmp=double(QstS[i])*double(potential[GrdPnt]);
		EnergySingular+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobNum[2];i++)
	{
		GrdPnt=IndexQmob[i];
		tmp=double(potential[GrdPnt])*double(Qmob[i])*double(dielectricChMob[i]);
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobDielBoarderNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarder[i];
		tmp=double(potential[GrdPnt])*double(QmobDielBoarder[i]);
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	
	

	totalEnergy=(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0);
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac = sqrt(fabs(oldSumSQ-SumSQ))/(fpoh*2.0);
	Dev=oldSumAbs-SumAbs;
	//if(oldSumSQ==0.0)fprintf(stdout,"PT Nit E dE RMSD dEabs\n");
	//if(verbose)fprintf(stdout,"PT %8d %17.12lg %17.12lg %17.12lg %17.12lg\n", iteration, totalEnergy, totalEnergy-oldTotEn,totalChange, Dev);
	ChargesEnergy[Counter]=totalEnergy;
	Counter++;
	if(Counter>=AvrOver)Counter=0;
	
	tmp=0.0;
	tmp2=0.0;
	for(i=0;i<AvrOver;i++)
	{
		tmp+=ChargesEnergy[i];
		tmp2+=ChargesEnergy[i]*ChargesEnergy[i];
	}
	tmp/=AvrOver;
	tmp2/=AvrOver;
	tmp2 = sqrt(fabs(tmp2-tmp*tmp));
	
	DbgPrint0("CSE_3       :E=%.16e Eq=%.10e Esing=%.10e EQmob=%.10e\n",(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));
	DbgPrint0("CSE_3       :E=%.16e Eq=%.10e\n",tmp,tmp2);
	
	oldSumSQ=SumSQ;
	oldSumAbs=SumAbs;
	oldTotEn=totalEnergy;
	World->SystemEnergy=totalEnergy;
	
	
	return EXIT_SUCCESS;
}
int compare_doubles( const void* a, const void* b ) 
{
	double* arg1 = (double*) a;
	double* arg2 = (double*) b;
	if( *arg1 < *arg2 ) return -1;
	else if( *arg1 == *arg2 ) return 0;
	else return 1;
}
int PoissonSolver::CalcSystemEnergy1(int iteration)
{
	int i,GrdPnt;
	int TotCharge=ChargeNum[2]+SingularNum[2]+QmobNum[2]+QmobDielBoarderNum[2],TotChargeCount;
	
	static double oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	double Dev=0.0;
	double tmp,SumSQ=0.0,SumAbs=0.0;
	double EnergyCharge = 0.0;
	double EnergySingular = 0.0;
	double EnergyQmob = 0.0;
	double OldTotalEnergy=totalEnergy;
	double  fpoh = 4.0*M_PI*World->GridScale;
	
	if(ChargesEnergy==NULL)
	{
		ChargesEnergy=new double[TotCharge];
	}
	TotChargeCount=0;
	for(i=0;i<ChargeNum[2];i++)
	{
		GrdPnt=IndexCharge[i];
		
		tmp=double(potential[GrdPnt])*double(Qst[i])*double(dielectricCh[i]);
		
		ChargesEnergy[TotChargeCount]=tmp;
		TotChargeCount++;
		
		EnergyCharge+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
		
	}
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		tmp=double(QstS[i])*double(potential[GrdPnt]);
		ChargesEnergy[TotChargeCount]=tmp;
		TotChargeCount++;
		
		EnergySingular+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobNum[2];i++)
	{
		GrdPnt=IndexQmob[i];
		tmp=double(potential[GrdPnt])*double(Qmob[i])*double(dielectricChMob[i]);
		ChargesEnergy[TotChargeCount]=tmp;
		TotChargeCount++;
		
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobDielBoarderNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarder[i];
		tmp=double(potential[GrdPnt])*double(QmobDielBoarder[i]);
		ChargesEnergy[TotChargeCount]=tmp;
		TotChargeCount++;
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	//sort ChargesEnergy
	qsort( ChargesEnergy, TotCharge, sizeof(double), compare_doubles );
	tmp=0.0;
	for(i=0;i<TotCharge;i++)
	{
		tmp+=ChargesEnergy[i];
	}
	DbgPrint0("CSE_1       :E=%.16e Eq=%.10e Esing=%.10e EQmob=%.10e\n",tmp/(fpoh*2.0),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));

	totalEnergy=(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0);
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac = sqrt(fabs(oldSumSQ-SumSQ))/(fpoh*2.0);
	Dev=oldSumAbs-SumAbs;
	//if(oldSumSQ==0.0)fprintf(stdout,"PT Nit E dE RMSD dEabs\n");
	//if(verbose)fprintf(stdout,"PT %8d %17.12lg %17.12lg %17.12lg %17.12lg\n", iteration, totalEnergy, totalEnergy-oldTotEn,totalChange, Dev);
	
	oldSumSQ=SumSQ;
	oldSumAbs=SumAbs;
	oldTotEn=totalEnergy;
	World->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
int PoissonSolver::CalcSystemEnergyAnalizer(int iteration)
{
	int i,GrdPnt;
	static long double oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	long double Dev=0.0;
	long double tmp,SumSQ=0.0,SumAbs=0.0;
	long double EnergyCharge = 0.0;
	long double EnergySingular = 0.0;
	long double EnergyQmob = 0.0;
	long double OldTotalEnergy=totalEnergy;
	long double ftmp, tmpAvr=0.0,tmpMin=0.0,tmpMax=0.0,tmpSD=0.0;
	long double  fpoh = 4.0*M_PI*World->GridScale;
	
	if(ChargeNum[2]>0)
	{
		tmp=fabs(potential[IndexCharge[0]]*Qst[0]*dielectricCh[0]);
		tmpMin=tmp;
		tmpMax=tmp;
	}
	
	for(i=0;i<ChargeNum[2];i++)
	{
		GrdPnt=IndexCharge[i];
		
		tmp=potential[GrdPnt]*Qst[i]*dielectricCh[i];
		EnergyCharge+=tmp;
		
		ftmp=fabs(tmp);
		tmpAvr+=ftmp;
		tmpSD+=ftmp*ftmp;
		
		if(ftmp<tmpMin)tmpMin=ftmp;
		if(ftmp>tmpMax)tmpMax=ftmp;
		
		//DbgPrint0("%.16Le %.16Le\n",tmp,ftmp);
		
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		tmp=QstS[i]*potential[GrdPnt];
		EnergySingular+=tmp;
		
		ftmp=fabs(tmp);
		tmpAvr+=ftmp;
		tmpSD+=ftmp*ftmp;
		
		if(ftmp<tmpMin)tmpMin=ftmp;
		if(ftmp>tmpMax)tmpMax=ftmp;
		
		//DbgPrint0("%.16Le %.16Le\n",tmp,ftmp);
		
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobNum[2];i++)
	{
		GrdPnt=IndexQmob[i];
		tmp=potential[GrdPnt]*Qmob[i]*dielectricChMob[i];
		
		ftmp=fabs(tmp);
		tmpAvr+=ftmp;
		tmpSD+=ftmp*ftmp;
		
		if(ftmp<tmpMin)tmpMin=ftmp;
		if(ftmp>tmpMax)tmpMax=ftmp;
		
		//DbgPrint0("%.16Le %.16Le\n",tmp,ftmp);
		
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobDielBoarderNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarder[i];
		tmp=potential[GrdPnt]*QmobDielBoarder[i];
		
		ftmp=fabs(tmp);
		tmpAvr+=ftmp;
		tmpSD+=ftmp*ftmp;
		
		if(ftmp<tmpMin)tmpMin=ftmp;
		if(ftmp>tmpMax)tmpMax=ftmp;
		
		//DbgPrint0("%.16Le %.16Le\n",tmp,ftmp);
		
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	tmpAvr=tmpAvr/(ChargeNum[2]+SingularNum[2]+QmobNum[2]+QmobDielBoarderNum[2]);
	tmpSD=tmpSD/(ChargeNum[2]+SingularNum[2]+QmobNum[2]+QmobDielBoarderNum[2]);
	tmpSD=sqrt(tmpSD-tmpAvr*tmpAvr);
	tmpAvr/=fpoh;
	tmpSD/=fpoh;
	tmpMin/=fpoh;
	tmpMax/=fpoh;
	DbgPrint0("fpoh %.7Le\n",fpoh);
	DbgPrint0("CSE_Analize :Avr(|Ei|)=%.5Le SD(|Ei|)=%.5Le Min(|Ei|)=%.10Le Max(|Ei|)=%.10Le\n",tmpAvr,tmpSD,tmpMin,tmpMax);
	DbgPrint0("CSE_LDouble :E=%.16Le Eq=%.10Le Esing=%.10Le EQmob=%.10lle\n",(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));

	totalEnergy=(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0);
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac = sqrt(fabs(oldSumSQ-SumSQ))/(fpoh*2.0);
	Dev=oldSumAbs-SumAbs;
	//if(oldSumSQ==0.0)fprintf(stdout,"PT Nit E dE RMSD dEabs\n");
	//if(verbose)fprintf(stdout,"PT %8d %17.12lg %17.12lg %17.12lg %17.12lg\n", iteration, totalEnergy, totalEnergy-oldTotEn,totalChange, Dev);
	
	oldSumSQ=SumSQ;
	oldSumAbs=SumAbs;
	oldTotEn=totalEnergy;
	World->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
int PoissonSolver::GuessNumberOfIteration()
{
	int GrdPnt;
	int i,j,k;
	int ix,iy,iz;
	float *sn[3];
	double tmp;
	for(i=0;i<3;i++)
	{
		sn[i]=new float [World->GridSize[i]];
		if(World->PBC[i])
		{
			for(j=0;j<World->GridSize[i];j++)
			{
				sn[i][j]=1.0/sqrt(float(World->GridSize[i]));
			}
		}
		else
		{
			sn[i][0]=0.0;
			sn[i][World->GridSize[i]-1]=0.0;
			for(j=1;j<World->GridSize[i]-1;j++)
			{
				tmp=M_PI*float(j)/float(World->GridSize[i]-1);
				sn[i][j]=sqrt(2.0)*sin(tmp)/sqrt(float(World->GridSize[i]-1));
			}
		}
	}
#ifndef PNPDOUBLE
	float *tmpPot=new float[GS_XYZ];
#else
	double *tmpPot=new double[GS_XYZ];
#endif
	
	//float *tmpPot2=new float[GS_XYZ];
	
	for(ix=0;ix<World->GridSize[0];ix++)
		for(iy=0;iy<World->GridSize[1];iy++)
			for(iz=0;iz<World->GridSize[2];iz++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		tmpPot[GrdPnt]=sn[0][ix]*sn[1][iy]*sn[2][iz];
		//tmpPot2[GrdPnt]=tmpPot[GrdPnt];
	}
#ifndef PNPDOUBLE
	float *PotPointer=potential;
#else
	double *PotPointer=potential;
#endif
	
	potential=tmpPot;
	int solverLoc;
	if(solver==Auto)
	{
		if(World->NIndexing==NULL)solverLoc=ArrayDirect;
		else solverLoc=NodeIndexBased;
	}
	else
	{
		solverLoc=solver;
	}
	
	if(solverLoc==NodeIndexBased)
	{
		float tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7;
		float om2 = 1.0;
		float om1 = 0.0;
		float om2d6 = 1.0/6.0;
		for(j=0;j<=1;j++){
			for(i=SingularNum[j];i<SingularNum[j+1];i++){
				GrdPnt=IndexSingular[i];
				tmp1 = potential[GrdPnt+1]*dielectricXS[i];
				tmp2 = potential[GrdPnt-1]*dielectricXmS[i];
				tmp3 = potential[GrdPnt+GS_X]*dielectricYS[i];
				tmp4 = potential[GrdPnt-GS_X]*dielectricYmS[i];
				tmp5 = potential[GrdPnt+GS_XY]*dielectricZS[i];
				tmp6 = potential[GrdPnt-GS_XY]*dielectricZmS[i];				
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3+tmp4;
				tmp3 = tmp5+tmp6;
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3;
				tmp1 = dielectricZSSUM[i]*(tmp1+tmp2);
				potential[GrdPnt] = potential[GrdPnt]*om1+tmp1;
			}
			for(i=ChargeNum[j];i<ChargeNum[j+1];i++){
				GrdPnt=IndexCharge[i];
				tmp1 = potential[GrdPnt+1]+potential[GrdPnt-1];
				tmp2 = potential[GrdPnt+GS_X]+potential[GrdPnt-GS_X];
				tmp3 = potential[GrdPnt+GS_XY]+potential[GrdPnt-GS_XY];
				tmp4 = potential[GrdPnt]*om1;
				tmp5 = om2d6*(tmp1+tmp2+tmp3);
				potential[GrdPnt] = tmp4+tmp5;
			}
			for(i=DielBoarderNum[j];i<DielBoarderNum[j+1];i++){
				GrdPnt=IndexDielBoarder[i];
				tmp1 = potential[GrdPnt+1]*dielectricXDB[i];
				tmp2 = potential[GrdPnt-1]*dielectricXmDB[i];
				tmp3 = potential[GrdPnt+GS_X]*dielectricYDB[i];
				tmp4 = potential[GrdPnt-GS_X]*dielectricYmDB[i];
				tmp5 = potential[GrdPnt+GS_XY]*dielectricZDB[i];
				tmp6 = potential[GrdPnt-GS_XY]*dielectricZmDB[i];
				tmp7 = potential[GrdPnt]*om1;
				potential[GrdPnt] = tmp7+dielectricZDBSUM[i]*(tmp1+tmp2+tmp3+tmp4+tmp5+tmp6);
			}
			for(i=NoSingularNum[j];i<NoSingularNum[j+1];i++){
				GrdPnt=IndexNoSingular[i];
				potential[GrdPnt] = potential[GrdPnt]*om1 + om2d6 * (potential[GrdPnt+1] + potential[GrdPnt-1] + potential[GrdPnt+GS_X] + potential[GrdPnt-GS_X] + potential[GrdPnt+GS_XY] + potential[GrdPnt-GS_XY]);
			}
			for(i=QmobNum[j];i<QmobNum[j+1];i++){
				GrdPnt=IndexQmob[i];
				tmp1 = potential[GrdPnt+1]+potential[GrdPnt-1];
				tmp2 = potential[GrdPnt+GS_X]+potential[GrdPnt-GS_X];
				tmp3 = potential[GrdPnt+GS_XY]+potential[GrdPnt-GS_XY];
				tmp4 = potential[GrdPnt]*om1;
				tmp5 = om2d6*(tmp1+tmp2+tmp3);
				potential[GrdPnt] = tmp4+tmp5;
			}
			for(i=QmobDielBoarderNum[j];i<QmobDielBoarderNum[j+1];i++)
			{
				GrdPnt=IndexQmobDielBoarder[i];
				tmp1 = potential[GrdPnt+1]*dielectricXQmobDB[i];
				tmp2 = potential[GrdPnt-1]*dielectricXmQmobDB[i];
				tmp3 = potential[GrdPnt+GS_X]*dielectricYQmobDB[i];
				tmp4 = potential[GrdPnt-GS_X]*dielectricYmQmobDB[i];
				tmp5 = potential[GrdPnt+GS_XY]*dielectricZQmobDB[i];
				tmp6 = potential[GrdPnt-GS_XY]*dielectricZmQmobDB[i];				
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3+tmp4;
				tmp3 = tmp5+tmp6;
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3;
				tmp1 = dielectricZQmobDBSUM[i]*(tmp1+tmp2);
				potential[GrdPnt] = potential[GrdPnt]*om1+tmp1;
			}
#	ifndef PNPDOUBLE
			World->BorderExchange(potential);
#	else
	
#	endif
			
		}
	}
	else if(solverLoc==ArrayDirect)
	{
		pnpError("Cann't estemate number of interaction for ArrayDirect solver\n");
		return 0;
	}
	else
	{
		pnpError("Unknown type of solver\n");
		return 0;
	}
	
	NodeIndex* NIndex=World->NIndexing->NIndex;
	unsigned int BlackAndWhiteMask=NodeIndexing::BlackAndWhiteMask;
	tmp=0.0;
	for(ix=0;ix<World->GridSize[0];ix++)
		for(iy=0;iy<World->GridSize[1];iy++)
			for(iz=0;iz<World->GridSize[2];iz++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		if(!(NIndex[GrdPnt]&BlackAndWhiteMask))
		{
		tmp+=potential[GrdPnt]*sn[0][ix]*sn[1][iy]*sn[2][iz];
		}
	}
	float spec=2.0*tmp;
	pnpPrint("Gauss-Seidel Spectral Radius: %f\n",spec);
	int iter=int(7.8/log(1.0 + sqrt(1.0-spec)));
	pnpPrint("Estimated Number of Iterations to Convergence: %d\n",iter);
	DeleteCVecArray(sn,3);
	delete [] tmpPot;
	potential=PotPointer;
	Relaxation=2.0/(1.0 + sqrt(1.0 - spec));
	pnpPrint("Estimated Relaxation Coefficient: %f\n",Relaxation);
	SetRelaxation(Relaxation);
	return iter;
}

////////////////////////////////////////////////////////////////////////////////
extern "C" int DoPoissonSolverOnCudaFloat(PoissonSolverOnCudaStruct* PS, PoissonSolverOnCudaParamStruct CudaParm);
PoissonSolverOnCuda::PoissonSolverOnCuda()
	: GenericSolver()
{
	InitZero();
}


PoissonSolverOnCuda::~PoissonSolverOnCuda()
{
	Clear();
}
int PoissonSolverOnCuda::InitZero()
{
	HaObject::SetName("PoissonSolverOnCuda");
	SolverStr.push_back("Auto");
	SolverStr.push_back("CPU");
	SolverStr.push_back("GPU");
	
	CudaParm.BS_X=96;
	CudaParm.BS_Y=2;
	CudaParm.BS_Z=1;
	CudaParm.BS_XY=CudaParm.BS_X*CudaParm.BS_Y;
	CudaParm.BS_XYZ=CudaParm.BS_XY*CudaParm.BS_Z;
	
	CudaParm.Qblock=32;
	CudaParm.DBblock=512;
	
	m_ContWorld=NULL;
	
	int i,j;
	PS.Relaxation=1.6;
	PS.MaxIterations=2000;
	PS.ConvergenceCheck=20;
	PS.Tolerance=0.0;
	PS.GS[0]=1;
	PS.GS[1]=1;
	PS.GS[2]=1;
	PS.spltGSWBC[0]=1;
	PS.spltGSWBC[1]=1;
	PS.spltGSWBC[2]=1;
	PS.GridScale=1.0;
	PS.spltGSWBC_X=1;
	PS.spltGSWBC_XY=1;
	PS.spltGSWBC_XYZ=1;
	PS.TotalEnergy=0.0;
	PS.AvrOverChecks=0;
	PS.TEavr=0.0;
	PS.stdevTE=0.0;
	
	for(i=0;i<8;i++)
	{
		PS.P[i]=NULL;
		PS.Qnum[i]=0;
		PS.Q[i]=NULL;
		PS.Qmult[i]=NULL;
		PS.Qpos[i]=NULL;
		
		PS.DielBordNum[i]=0;
		PS.DielBordPos[i]=NULL;
		for(j=0;j<6;j++)
			PS.DielMult[i][j]=NULL;
	}
	return EXIT_SUCCESS;
}
int PoissonSolverOnCuda::Clear()
{
	int i,j;
	for(i=0;i<8;i++)
	{
		DeleteCArray(PS.P[i]);
		DeleteCArray(PS.Q[i]);
		DeleteCArray(PS.Qmult[i]);
		DeleteCArray(PS.Qpos[i]);
		if(PS.DielBordPos[i]!=NULL)
		{
			delete [] PS.DielBordPos[i];
			PS.DielBordPos[i]=NULL;
		}
		//DeleteCArray(PS.DielBordPos[i]);
		for(j=0;j<6;j++)
			DeleteCArray(PS.DielMult[i][j]);
	}
	return EXIT_SUCCESS;
}
int PoissonSolverOnCuda::SaveXML(TiXmlElement* Elt, HaContext* p_ctxt )
{
	return EXIT_SUCCESS;
}
int PoissonSolverOnCuda::LoadXML(const TiXmlElement* Elt, HaContext* p_ctxt )
{
	if(strncmp(HaObject::GetCStrName(),Elt->Value(),13))
	{
		fprintf(stderr,"ERROR: Wrong XML Element %s, expecting %s\n",Elt->Value(),HaObject::GetCStrName());
		return EXIT_FAILURE;
	}
	//Read Primary Parameters
	Elt->GetIntAttribute("MaxIterations",&PS.MaxIterations);
	Elt->GetFloatAttribute("Relaxation",&PS.Relaxation);
	Elt->GetIntAttribute("ConvergenceCheck",&PS.ConvergenceCheck);
	
	if(Elt->GetStdStrIndex("Solver",&solver,SolverStr)!=EXIT_SUCCESS)solver=2;
	if(solver==0)solver=2;
	return EXIT_SUCCESS;
}
int PoissonSolverOnCuda::SetContWorld(ContWorld* _world)
{
	m_ContWorld=_world;
	PS.GS[0]=m_ContWorld->GridSizeGlobal[0]-2;
	PS.GS[1]=m_ContWorld->GridSizeGlobal[1]-2;
	PS.GS[2]=m_ContWorld->GridSizeGlobal[2]-2;
	PS.spltGSWBC[0]=PS.GS[0]/2+CUDAXTRAX;
	PS.spltGSWBC[1]=PS.GS[1]/2+2;
	PS.spltGSWBC[2]=PS.GS[2]/2+2;
	PS.spltGSWBC_X=PS.spltGSWBC[0];
	PS.spltGSWBC_XY=PS.spltGSWBC[0]*PS.spltGSWBC[1];
	PS.spltGSWBC_XYZ=PS.spltGSWBC[0]*PS.spltGSWBC[1]*PS.spltGSWBC[2];
	PS.GridScale=m_ContWorld->GridScale;
	return EXIT_SUCCESS;
}
int PoissonSolverOnCuda::ShowParameters()
{
	return EXIT_SUCCESS;
}
int PoissonSolverOnCuda::ShowProperties()
{
	return EXIT_SUCCESS;
}
int PoissonSolverOnCuda::GuessNumberOfIteration()
{
	pnpPrint("GuessNumberOfIteration\n");
	int GrdPnt;
	int i,j,k;
	int ix,iy,iz;
	float *sn[3];
	double tmp;
	for(i=0;i<3;i++)
	{
		sn[i]=new float [m_ContWorld->GridSize[i]];
		if(m_ContWorld->PBC[i])
		{
			for(j=0;j<m_ContWorld->GridSize[i];j++)
			{
				sn[i][j]=1.0/sqrt(float(m_ContWorld->GridSize[i]));
			}
		}
		else
		{
			sn[i][0]=0.0;
			sn[i][m_ContWorld->GridSize[i]-1]=0.0;
			for(j=1;j<m_ContWorld->GridSize[i]-1;j++)
			{
				tmp=M_PI*float(j)/float(m_ContWorld->GridSize[i]-1);
				sn[i][j]=sqrt(2.0)*sin(tmp)/sqrt(float(m_ContWorld->GridSize[i]-1));
			}
		}
	}
	
	//set split potential
	int spltGSWBC_XYZ=PS.spltGSWBC[0]*PS.spltGSWBC[1]*PS.spltGSWBC[2];
	for(i=0;i<8;i++)
	{
		for(j=0;j<spltGSWBC_XYZ;j++)
			PS.P[i][j]=0.0f;
	}
	int tx,ty,tz;
	int PotGS_X=m_ContWorld->GS_X;
	int PotGS_XY=m_ContWorld->GS_XY;
	//set split potential
	for(tz=1;tz<PS.spltGSWBC[2]-1;tz++)
		for(ty=1;ty<PS.spltGSWBC[1]-1;ty++)
			for(tx=0;tx<PS.spltGSWBC[0]-CUDAXTRAX;tx++)
	{
		ix = 2*tx+1;
		iy = 2*ty-1;
		iz = 2*tz-1;
		
		i=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
		//GrdPnt0=ix+iy*PotGS_X+iz*PotGS_XY;
		PS.P[0][i]=sn[0][ix]*sn[1][iy]*sn[2][iz];
		PS.P[1][i]=sn[0][ix+1]*sn[1][iy]*sn[2][iz];
		PS.P[2][i]=sn[0][ix]*sn[1][iy+1]*sn[2][iz];
		PS.P[3][i]=sn[0][ix+1]*sn[1][iy+1]*sn[2][iz];
		PS.P[4][i]=sn[0][ix]*sn[1][iy]*sn[2][iz+1];
		PS.P[5][i]=sn[0][ix+1]*sn[1][iy]*sn[2][iz+1];
		PS.P[6][i]=sn[0][ix]*sn[1][iy+1]*sn[2][iz+1];
		PS.P[7][i]=sn[0][ix+1]*sn[1][iy+1]*sn[2][iz+1];
	}
	
	//solve once
	int iter=PS.MaxIterations;
	PS.MaxIterations=1;
	int locQnum[8];
	for(i=0;i<8;i++)
	{
		locQnum[i]=PS.Qnum[i];
		PS.Qnum[i]=0;
	}
	DoPoissonSolverCudaOnCPU();
	for(i=0;i<8;i++)
	{
		PS.Qnum[i]=locQnum[i];
	}
	PS.MaxIterations=iter;
	
	//calculating relaxation and number of iterations
	tmp=0.0;
	
	for(tz=1;tz<PS.spltGSWBC[2]-1;tz++)
		for(ty=1;ty<PS.spltGSWBC[1]-1;ty++)
			for(tx=0;tx<PS.spltGSWBC[0]-CUDAXTRAX;tx++)
	{
		ix = 2*tx+1;
		iy = 2*ty-1;
		iz = 2*tz-1;
		
		i=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
		//GrdPnt0=ix+iy*PotGS_X+iz*PotGS_XY;
// 		tmp+=PS.P[0][i]*sn[0][ix]*sn[1][iy]*sn[2][iz];
// 		tmp+=PS.P[3][i]*sn[0][ix+1]*sn[1][iy+1]*sn[2][iz];
// 		tmp+=PS.P[5][i]*sn[0][ix+1]*sn[1][iy]*sn[2][iz+1];
// 		tmp+PS.P[6][i]*sn[0][ix]*sn[1][iy+1]*sn[2][iz+1];
		tmp+=PS.P[1][i]*sn[0][ix+1]*sn[1][iy]*sn[2][iz];
		tmp+=PS.P[2][i]*sn[0][ix]*sn[1][iy+1]*sn[2][iz];
		tmp+=PS.P[4][i]*sn[0][ix]*sn[1][iy]*sn[2][iz+1];
		tmp+=PS.P[7][i]*sn[0][ix+1]*sn[1][iy+1]*sn[2][iz+1];
		
	}
	float spec=2.0*tmp;
	pnpPrint("Gauss-Seidel Spectral Radius: %f\n",spec);
	iter=int(7.8/log(1.0 + sqrt(1.0-spec)));
	
	DeleteCVecArray(sn,3);
	PS.Relaxation=2.0/(1.0 + sqrt(1.0 - spec));
	for(i=0;i<8;i++)
	{
		for(j=0;j<spltGSWBC_XYZ;j++)
			PS.P[i][j]=0.0f;
	}
	pnpPrint("Estimated Number of Iterations to Convergence: %d\n",iter);
	pnpPrint("Estimated Relaxation Coefficient: %f\n",PS.Relaxation);
	
	//iter=10000;
	//PS.Relaxation=1.9;
	//pnpPrint("Estimated Number of Iterations to Convergence: %d\n",iter);
	//pnpPrint("Estimated Relaxation Coefficient: %f\n",PS.Relaxation);
	
	if(PS.MaxIterations<0)
	{
		PS.MaxIterations=iter;
		if(PS.MaxIterations<60)
		{
			pnpPrint("Estimated Number of Iterations <60 will set to 60\n");
			PS.MaxIterations=60;
		}
	}
	return 1;
}
int PoissonSolverOnCuda::InitSolver()
{
	DbgPrint0("PoissonSolverOnCuda::InitSolver()\n");
	
	bool bGuessNumberOfIteration=false;
	if(PS.Relaxation<0.0)
	{
		bGuessNumberOfIteration=true;
		PS.Relaxation=1.0;
	}
	
	int i,j,k,GrdPnt;
	int ix,iy,iz;
	int t,tx,ty,tz;
	float ftmp1,ftmp2,ftmp3;
	float om2=PS.Relaxation;
	float om2d6=PS.Relaxation/6.0;
	NodeIndexing* NIndexing=m_ContWorld->NIndexing;
	NodeIndex* NIndex=NIndexing->NIndex;
	float *Eps=NIndexing->Eps;
	
	NodeIndex specChargeMask=NodeIndexing::ChargeMask;
	NodeIndex DielBoarderMask=NodeIndexing::DielBoarderMask;
	
	int spltGSWBC_XYZ=PS.spltGSWBC[0]*PS.spltGSWBC[1]*PS.spltGSWBC[2];
	//allocate host split potential
	for(i=0;i<8;i++)
	{
		PS.P[i]=new float[spltGSWBC_XYZ];
		for(j=0;j<spltGSWBC_XYZ;j++)
			PS.P[i][j]=0.0f;
	}
	
	
	int pGS_X=m_ContWorld->GS_X;
	int pGS_Y=m_ContWorld->GS_Y;
	int pGS_Z=m_ContWorld->GS_Z;
	int pGS_XY=pGS_X*pGS_Y;
	int pGS_XYZ=pGS_X*pGS_Y*pGS_Z;
	
	int izgrid,iygrid;
	//cound Q and BielBorder
	for(i=0;i<8;i++)
	{
		PS.Qnum[i]=0;
		PS.DielBordNum[i]=0;
	}
	int QCount=0;
	int DBCount=0;
	for(iz=1;iz<pGS_Z-1;iz++)
	{
		izgrid = iz*pGS_XY;
		for(iy=1;iy<pGS_Y-1;iy++)
		{
			iygrid = izgrid+iy*pGS_X;
			for(ix=1;ix<pGS_X-1;ix++)
			{
				GrdPnt = iygrid+ix;
				if(NIndex[GrdPnt]&specChargeMask)
				{
					i=((ix - 1)%2)+((iy + 1)%2)*2+((iz + 1)%2)*4;
					PS.Qnum[i]++;
					QCount++;
				}
				if(NIndex[GrdPnt]&DielBoarderMask)
				{
					i=((ix - 1)%2)+((iy + 1)%2)*2+((iz + 1)%2)*4;
					PS.DielBordNum[i]++;
					DBCount++;
				}
			}
		}
	}
	for(i=0;i<8;i++)
		DbgPrint0("PS.Qnum[%d]=%d\n",i,PS.Qnum[i]);
	DbgPrint0("QCount %d\n",QCount);
	for(i=0;i<8;i++)
		DbgPrint0("PS.DielBordNum[%d]=%d\n",i,PS.DielBordNum[i]);
	DbgPrint0("DBCount %d\n",DBCount);
	//allocate split charge and diel border
	int iQcount[8];
	int iDBcount[8];
	for(i=0;i<8;i++)
	{
		if(PS.Qnum[i]>0)
		{
			PS.Q[i]=new float[PS.Qnum[i]];
			PS.Qpos[i]=new int[PS.Qnum[i]];
			PS.Qmult[i]=new float[PS.Qnum[i]];
			for(j=0;j<PS.Qnum[i];j++)
			{
				PS.Q[i][j]=0.0f;
				PS.Qpos[i][j]=-1;
				PS.Qmult[i][j]=0.0f;
			}
		}
		else
		{
			PS.Q[i]=NULL;
			PS.Qpos[i]=NULL;
			PS.Qmult[i]=NULL;
		}
		if(PS.DielBordNum[i]>0)
		{
			PS.DielBordPos[i]=new int[PS.DielBordNum[i]];
			for(j=0;j<PS.DielBordNum[i];j++)
				PS.DielBordPos[i][j]=-1;
			for(k=0;k<6;k++)
			{
				PS.DielMult[i][k]=new float[PS.DielBordNum[i]];
				for(j=0;j<PS.DielBordNum[i];j++)
					PS.DielMult[i][k][j]=0.0;
			}
		}
		else
		{
			PS.DielBordPos[i]==NULL;
			for(k=0;k<6;k++)
				PS.DielMult[i][k]=NULL;
		}
		iQcount[i]=0;
		iDBcount[i]=0;
	}
	//fill split charge and diel border
	//int count=0;
	QCount=0;
	for(iz=1;iz<pGS_Z-1;iz++)
	{
		izgrid = iz*pGS_XY;
		for(iy=1;iy<pGS_Y-1;iy++)
		{
			iygrid = izgrid+iy*pGS_X;
			for(ix=1;ix<pGS_X-1;ix++)
			{
				GrdPnt = iygrid+ix;
				if(NIndex[GrdPnt]&specChargeMask)
				{
					tx=(ix - 1)/2;
					ty=(iy + 1)/2;
					tz=(iz + 1)/2;
					t=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
					i=((ix - 1)%2)+((iy + 1)%2)*2+((iz + 1)%2)*4;
					
					PS.Qpos[i][iQcount[i]]=t;
					if(NIndex[GrdPnt]&DielBoarderMask)
					{
						ftmp1=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft] + Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft] + Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft] + Eps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft] + Eps[(NIndex[GrdPnt-pGS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft] + Eps[(NIndex[GrdPnt-pGS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
						PS.Qmult[i][iQcount[i]]=om2/ftmp1;
						PS.Q[i][iQcount[i]]=PS.Qmult[i][iQcount[i]]*NIndexing->Q[QCount];
						
					}
					else
					{
						PS.Qmult[i][iQcount[i]]=om2d6/Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
						PS.Q[i][iQcount[i]]=PS.Qmult[i][iQcount[i]]*NIndexing->Q[QCount];
					}
					
					
					iQcount[i]++;
					QCount++;
				}
				if(NIndex[GrdPnt]&DielBoarderMask)
				{
					tx=(ix - 1)/2;
					ty=(iy + 1)/2;
					tz=(iz + 1)/2;
					t=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
					i=((ix - 1)%2)+((iy + 1)%2)*2+((iz + 1)%2)*4;
					
					PS.DielBordPos[i][iDBcount[i]]=t;
					
					PS.DielMult[i][PlusX][iDBcount[i]] =Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
					PS.DielMult[i][PlusY][iDBcount[i]] =Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
					PS.DielMult[i][PlusZ][iDBcount[i]] =Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
					PS.DielMult[i][MinusX][iDBcount[i]]=Eps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
					PS.DielMult[i][MinusY][iDBcount[i]]=Eps[(NIndex[GrdPnt-pGS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
					PS.DielMult[i][MinusZ][iDBcount[i]]=Eps[(NIndex[GrdPnt-pGS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
					
					ftmp1=PS.DielMult[i][PlusX][iDBcount[i]]+PS.DielMult[i][MinusX][iDBcount[i]]+PS.DielMult[i][PlusY][iDBcount[i]]+PS.DielMult[i][MinusY][iDBcount[i]]+PS.DielMult[i][PlusZ][iDBcount[i]]+PS.DielMult[i][MinusZ][iDBcount[i]];
					ftmp1=PS.Relaxation/ftmp1;
					
					for(k=0;k<6;k++)
					{
						PS.DielMult[i][k][iDBcount[i]] = ftmp1*PS.DielMult[i][k][iDBcount[i]]-om2d6;
					}
// 					if(count<10)
// 					{
// 						printf("[%f %f %f %f %f %f]  %f %f\n",PS.DielMult[i][PlusX][iDBcount[i]], PS.DielMult[i][MinusX][iDBcount[i]], PS.DielMult[i][PlusY][iDBcount[i]], PS.DielMult[i][MinusY][iDBcount[i]], PS.DielMult[i][PlusZ][iDBcount[i]], PS.DielMult[i][MinusZ][iDBcount[i]], om2d6, ftmp1);
// 						count++;
// 					}
					iDBcount[i]++;
				}
			}
		}
	}
	if(bGuessNumberOfIteration)
	{
		GuessNumberOfIteration();
		InitSolverSetIntArrRel();
	}
	return EXIT_SUCCESS;
}
int PoissonSolverOnCuda::InitSolverSetIntArrRel()
{
	int i,j,k,GrdPnt;
	int ix,iy,iz;
	int t,tx,ty,tz;
	float ftmp1,ftmp2,ftmp3;
	float om2=PS.Relaxation;
	float om2d6=PS.Relaxation/6.0;
	NodeIndexing* NIndexing=m_ContWorld->NIndexing;
	NodeIndex* NIndex=NIndexing->NIndex;
	float *Eps=NIndexing->Eps;
	
	NodeIndex specChargeMask=NodeIndexing::ChargeMask;
	NodeIndex DielBoarderMask=NodeIndexing::DielBoarderMask;
	
	int spltGSWBC_XYZ=PS.spltGSWBC[0]*PS.spltGSWBC[1]*PS.spltGSWBC[2];
	
	int pGS_X=m_ContWorld->GS_X;
	int pGS_Y=m_ContWorld->GS_Y;
	int pGS_Z=m_ContWorld->GS_Z;
	int pGS_XY=pGS_X*pGS_Y;
	int pGS_XYZ=pGS_X*pGS_Y*pGS_Z;
	
	int izgrid,iygrid;
	
	int QCount=0;
	int DBCount=0;
	
	int iQcount[8];
	int iDBcount[8];
	for(i=0;i<8;i++)
	{
		iQcount[i]=0;
		iDBcount[i]=0;
	}
	//fill split charge and diel border
	QCount=0;
	for(iz=1;iz<pGS_Z-1;iz++)
	{
		izgrid = iz*pGS_XY;
		for(iy=1;iy<pGS_Y-1;iy++)
		{
			iygrid = izgrid+iy*pGS_X;
			for(ix=1;ix<pGS_X-1;ix++)
			{
				GrdPnt = iygrid+ix;
				if(NIndex[GrdPnt]&specChargeMask)
				{
					tx=(ix - 1)/2;
					ty=(iy + 1)/2;
					tz=(iz + 1)/2;
					t=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
					i=((ix - 1)%2)+((iy + 1)%2)*2+((iz + 1)%2)*4;
					
					PS.Qpos[i][iQcount[i]]=t;
					if(NIndex[GrdPnt]&DielBoarderMask)
					{
						ftmp1=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft] + Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft] + Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft] + Eps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft] + Eps[(NIndex[GrdPnt-pGS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft] + Eps[(NIndex[GrdPnt-pGS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
						PS.Qmult[i][iQcount[i]]=om2/ftmp1;
						PS.Q[i][iQcount[i]]=PS.Qmult[i][iQcount[i]]*NIndexing->Q[QCount];
					}
					else
					{
						PS.Qmult[i][iQcount[i]]=om2d6/Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
						PS.Q[i][iQcount[i]]=PS.Qmult[i][iQcount[i]]*NIndexing->Q[QCount];
					}
					
					
					iQcount[i]++;
					QCount++;
				}
				if(NIndex[GrdPnt]&DielBoarderMask)
				{
					tx=(ix - 1)/2;
					ty=(iy + 1)/2;
					tz=(iz + 1)/2;
					t=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
					i=((ix - 1)%2)+((iy + 1)%2)*2+((iz + 1)%2)*4;
					
					PS.DielBordPos[i][iDBcount[i]]=t;
					
					PS.DielMult[i][PlusX][iDBcount[i]] =Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
					PS.DielMult[i][PlusY][iDBcount[i]] =Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
					PS.DielMult[i][PlusZ][iDBcount[i]] =Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
					PS.DielMult[i][MinusX][iDBcount[i]]=Eps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
					PS.DielMult[i][MinusY][iDBcount[i]]=Eps[(NIndex[GrdPnt-pGS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
					PS.DielMult[i][MinusZ][iDBcount[i]]=Eps[(NIndex[GrdPnt-pGS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
					
					ftmp1=PS.DielMult[i][PlusX][iDBcount[i]]+PS.DielMult[i][MinusX][iDBcount[i]]+PS.DielMult[i][PlusY][iDBcount[i]]+PS.DielMult[i][MinusY][iDBcount[i]]+PS.DielMult[i][PlusZ][iDBcount[i]]+PS.DielMult[i][MinusZ][iDBcount[i]];
					ftmp1=PS.Relaxation/ftmp1;
					
					for(k=0;k<6;k++)
					{
						PS.DielMult[i][k][iDBcount[i]] = ftmp1*PS.DielMult[i][k][iDBcount[i]]-om2d6;
					}
					iDBcount[i]++;
				}
			}
		}
	}
	return EXIT_SUCCESS;
}
int PoissonSolverOnCuda::Solve()
{
	DbgPrint0("PoissonSolverOnCuda::Solve()\n");
	
	if(m_ContWorld->Potential == NULL)
	{
		m_ContWorld->Potential = new float[m_ContWorld->GS_XYZ];
		for(int i=0;i<m_ContWorld->GS_XYZ;i++)
			m_ContWorld->Potential[i]= 0.0;
	}
	
	CopyPotBC2SplitPot();
	CopyIntPot2SplitPot();
	if(solver==1)
		DoPoissonSolverCudaOnCPU();
	else
	{
		int halfGS_X=(m_ContWorld->GS_X-2)/2;
		int FindBS=0;
		if(CudaParm.BS_X==0)
		{
			FindBS=1;
		}
		else
		{
			if(halfGS_X%CudaParm.BS_X!=0)
				FindBS=1;
		}
		if(FindBS)
		{
			for(CudaParm.BS_X=128;CudaParm.BS_X>=32;CudaParm.BS_X=CudaParm.BS_X-32)
			{
				DbgPrint0("halfGS_XCudaParm.BS_X %d %d %d %d\n",m_ContWorld->GS_X,halfGS_X,CudaParm.BS_X,halfGS_X%CudaParm.BS_X);
				if(halfGS_X%CudaParm.BS_X==0)
				{
					DbgPrint0("BS_X=%d and it is optimal\n",CudaParm.BS_X);
					break;
				}
			}
			if(halfGS_X%32!=0)
			{
				if(halfGS_X<96)
				{
					CudaParm.BS_X=halfGS_X;
					DbgPrint0("BS_X=%d and it is NOT optimal\n",CudaParm.BS_X);
				}
				else
				{
					pnpError("Can not find block size for GPU, change GridSize of the system\n");
				}
			}
			
			CudaParm.BS_Y=2;
			CudaParm.BS_Z=1;
		}
		CudaParm.BS_XY=CudaParm.BS_X*CudaParm.BS_Y;
		CudaParm.BS_XYZ=CudaParm.BS_XY*CudaParm.BS_Z;
#if defined(WITH_CUDA)
		DoPoissonSolverOnCudaFloat(&PS,CudaParm);
#endif
	}
	CopySplitPot2IntPot();
	
	
// 	CopyPotBC2SplitPot();
// 	int Qblocks[6]={16,32,64,128,256,512};
// 	for(int i=0;i<6;i++)
// 	{
// 		CudaParm.Qblock=Qblocks[i];
// 		CopyIntPot2SplitPot();
// 		DoPoissonSolverOnCuda(PS,CudaParm);
// 	}
	//calculate energy
	CalcSystemEnergy();
	pnpPrintGroup0("<PoissonFinal E=\"%.10e\"/>\n", m_ContWorld->SystemEnergy);
	return EXIT_SUCCESS;
}
int PoissonSolverOnCuda::CalcSystemEnergy()
{
	int i,j;
	static double oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	double Dev=0.0;
	double tmp,SumSQ=0.0,SumAbs=0.0;
	double EnergyCharge = 0.0;
	double EnergySingular = 0.0;
	double EnergyQmob = 0.0;
	double totalEnergy=0.0;
	double fpoh = 4.0*M_PI*m_ContWorld->GridScale;
	float dphi=0.0,maxdphi=0.0;
	int GrdPntMaxDPhi=-1;
	float om2d6=PS.Relaxation/6.0;
	for(i=0;i<8;i++)
	{
		for(j=0;j<PS.Qnum[i];j++)
		{
			tmp=double(PS.P[i][PS.Qpos[i][j]])*double(PS.Q[i][j])/double(PS.Qmult[i][j]);
			EnergyCharge+=tmp;
		}
	}
	totalEnergy=EnergyCharge/(fpoh*2.0);
	m_ContWorld->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
int PoissonSolverOnCuda::CopyPotBC2SplitPot()
{
	float* Pot=m_ContWorld->Potential;
	int GrdPnt0,ix,iy,iz;
	int i,tx,ty,tz;
	int PotGS_X=m_ContWorld->GS_X;
	int PotGS_XY=m_ContWorld->GS_XY;
	//copy BC Potential to split potential
	//XY
	for(ty=1;ty<PS.spltGSWBC[1]-1;ty++)
		for(tx=0;tx<PS.spltGSWBC[0]-CUDAXTRAX;tx++)
	{
		tz=0;
		ix = 2*tx+1;iy = 2*ty-1;iz = 2*tz-1;
		i=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
		GrdPnt0=ix+iy*PotGS_X+iz*PotGS_XY;
		
		PS.P[4][i]=Pot[GrdPnt0+PotGS_XY];
		PS.P[5][i]=Pot[GrdPnt0+PotGS_XY+1];
		PS.P[6][i]=Pot[GrdPnt0+PotGS_XY+PotGS_X];
		PS.P[7][i]=Pot[GrdPnt0+PotGS_XY+1+PotGS_X];
			
		tz=PS.spltGSWBC[2]-1;
		ix = 2*tx+1;iy = 2*ty-1;iz = 2*tz-1;	
		i=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
		GrdPnt0=ix+iy*PotGS_X+iz*PotGS_XY;
		
		PS.P[0][i]=Pot[GrdPnt0];
		PS.P[1][i]=Pot[GrdPnt0+1];
		PS.P[2][i]=Pot[GrdPnt0+PotGS_X];
		PS.P[3][i]=Pot[GrdPnt0+1+PotGS_X];
		
	}
	//ZX
	for(tz=1;tz<PS.spltGSWBC[2]-1;tz++)
		for(tx=0;tx<PS.spltGSWBC[0]-CUDAXTRAX;tx++)
	{
		ty=0;
		ix = 2*tx+1;iy = 2*ty-1;iz = 2*tz-1;
		i=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
		GrdPnt0=ix+iy*PotGS_X+iz*PotGS_XY;
			
		PS.P[2][i]=Pot[GrdPnt0+PotGS_X];
		PS.P[3][i]=Pot[GrdPnt0+1+PotGS_X];
		PS.P[6][i]=Pot[GrdPnt0+PotGS_XY+PotGS_X];
		PS.P[7][i]=Pot[GrdPnt0+PotGS_XY+1+PotGS_X];
			
		ty=PS.spltGSWBC[1]-1;
		ix = 2*tx+1;iy = 2*ty-1;iz = 2*tz-1;
		i=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
		GrdPnt0=ix+iy*PotGS_X+iz*PotGS_XY;
			
		PS.P[0][i]=Pot[GrdPnt0];
		PS.P[1][i]=Pot[GrdPnt0+1];
		PS.P[4][i]=Pot[GrdPnt0+PotGS_XY];
		PS.P[5][i]=Pot[GrdPnt0+PotGS_XY+1];
	}
	//YZ
	for(tz=1;tz<PS.spltGSWBC[2]-1;tz++)
		for(ty=1;ty<PS.spltGSWBC[1]-1;ty++)
	{
		tx=-1;
		ix = 2*tx+1;iy = 2*ty-1;iz = 2*tz-1;
		i=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
		GrdPnt0=ix+iy*PotGS_X+iz*PotGS_XY;
		PS.P[1][i]=Pot[GrdPnt0+1];
		PS.P[3][i]=Pot[GrdPnt0+1+PotGS_X];
		PS.P[5][i]=Pot[GrdPnt0+PotGS_XY+1];
		PS.P[7][i]=Pot[GrdPnt0+PotGS_XY+1+PotGS_X];
			
		tx=PS.spltGSWBC[0]-CUDAXTRAX;
		ix = 2*tx+1;iy = 2*ty-1;iz = 2*tz-1;
		i=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
		GrdPnt0=ix+iy*PotGS_X+iz*PotGS_XY;
		PS.P[0][i]=Pot[GrdPnt0];
		PS.P[2][i]=Pot[GrdPnt0+PotGS_X];
		PS.P[4][i]=Pot[GrdPnt0+PotGS_XY];
		PS.P[6][i]=Pot[GrdPnt0+PotGS_XY+PotGS_X];
			
	}
	return EXIT_SUCCESS;
}
int PoissonSolverOnCuda::CopyIntPot2SplitPot()
{
	float* Pot=m_ContWorld->Potential;
	int GrdPnt0,ix,iy,iz;
	int i,tx,ty,tz;
	int PotGS_X=m_ContWorld->GS_X;
	int PotGS_XY=m_ContWorld->GS_XY;
	//copy internal Potential to split potential
	for(tz=1;tz<PS.spltGSWBC[2]-1;tz++)
		for(ty=1;ty<PS.spltGSWBC[1]-1;ty++)
			for(tx=0;tx<PS.spltGSWBC[0]-CUDAXTRAX;tx++)
	{
		ix = 2*tx+1;
		iy = 2*ty-1;
		iz = 2*tz-1;
			
		i=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
		GrdPnt0=ix+iy*PotGS_X+iz*PotGS_XY;
		PS.P[0][i]=Pot[GrdPnt0];
		PS.P[1][i]=Pot[GrdPnt0+1];
		PS.P[2][i]=Pot[GrdPnt0+PotGS_X];
		PS.P[3][i]=Pot[GrdPnt0+1+PotGS_X];
		PS.P[4][i]=Pot[GrdPnt0+PotGS_XY];
		PS.P[5][i]=Pot[GrdPnt0+PotGS_XY+1];
		PS.P[6][i]=Pot[GrdPnt0+PotGS_XY+PotGS_X];
		PS.P[7][i]=Pot[GrdPnt0+PotGS_XY+1+PotGS_X];
	}
	return EXIT_SUCCESS;
}
int PoissonSolverOnCuda::CopySplitPot2IntPot()
{
	float* Pot=m_ContWorld->Potential;
	int GrdPnt0,ix,iy,iz;
	int i,tx,ty,tz;
	int PotGS_X=m_ContWorld->GS_X;
	int PotGS_XY=m_ContWorld->GS_XY;
	//copy internal part of split potential to internal part of Potential
	for(tz=1;tz<PS.spltGSWBC[2]-1;tz++)
		for(ty=1;ty<PS.spltGSWBC[1]-1;ty++)
			for(tx=0;tx<PS.spltGSWBC[0]-CUDAXTRAX;tx++)
	{
		ix = 2*tx+1;
		iy = 2*ty-1;
		iz = 2*tz-1;
			
		i=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
		GrdPnt0=ix+iy*PotGS_X+iz*PotGS_XY;
		Pot[GrdPnt0]=PS.P[0][i];
		Pot[GrdPnt0+1]=PS.P[1][i];
		Pot[GrdPnt0+PotGS_X]=PS.P[2][i];
		Pot[GrdPnt0+1+PotGS_X]=PS.P[3][i];
		Pot[GrdPnt0+PotGS_XY]=PS.P[4][i];
		Pot[GrdPnt0+PotGS_XY+1]=PS.P[5][i];
		Pot[GrdPnt0+PotGS_XY+PotGS_X]=PS.P[6][i];
		Pot[GrdPnt0+PotGS_XY+1+PotGS_X]=PS.P[7][i];
	}
	return EXIT_SUCCESS;
}

int PoissonSolverOnCuda::DoPoissonSolverCudaOnCPU()
{
	pnpPrintGroup0("Solving Poisson Equation using algorithm made for CUDA\n");
	int i,j,k,t,tx,ty,tz;
	int pitchX=PS.spltGSWBC_X;
	int pitchXY=PS.spltGSWBC_XY;
	int spltGSWBC_X=PS.spltGSWBC[0];
	int spltGSWBC_Y=PS.spltGSWBC[1];
	int spltGSWBC_Z=PS.spltGSWBC[2];
	int MaxIterations=PS.MaxIterations;
	int ConvergenceCheck=PS.ConvergenceCheck;
	
	
	int iteration;
	float om2 = PS.Relaxation;
	float om1 = 1.0-om2;
	float om2d6 = om2/6.0;
	
	float *d_P0=PS.P[0];
	float *d_P1=PS.P[1];
	float *d_P2=PS.P[2];
	float *d_P3=PS.P[3];
	float *d_P4=PS.P[4];
	float *d_P5=PS.P[5];
	float *d_P6=PS.P[6];
	float *d_P7=PS.P[7];
	
	float ftmp1,ftmp2;
	double dtmp1,dtmp2;
	
	double fpoh = 4.0*M_PI*m_ContWorld->GridScale;
	float ftmp;
	double totalEnergy=0.0;
	double totalChange=0.0;
	double relativeChange=0.0;
	double ConvFac=0.0;
	double EnergyCharge = 0.0;
	
	int d_Qnum[8];
	float *d_Q[8];
	int *d_Qpos[8];
	
	float* dc_DielMult[48];
	int dc_DBnum[8];
	int* dc_DBpos[8];
	for(i=0;i<8;i++)
	{
		dc_DBpos[i]=NULL;
		for(k=0;k<6;k++)
			dc_DielMult[i*6+k]=NULL;
		if(PS.DielBordNum[i]>0)
		{
			for(k=0;k<6;k++)
			{
				dc_DielMult[i*6+k]=PS.DielMult[i][k];
			}
			dc_DBpos[i]=PS.DielBordPos[i];
		}
		dc_DBnum[i]=PS.DielBordNum[i];
	}
	
			
	
	for(i=0;i<8;i++)
	{
		d_Qnum[i]=PS.Qnum[i];
		d_Q[i]=PS.Q[i];
		d_Qpos[i]=PS.Qpos[i];
	}
	
	
	for(i=0;i<8;i++)
		DbgPrint0("PS.Qnum[%d]=%d\n",i,d_Qnum[i]);
	for(iteration=1;iteration<=MaxIterations;iteration++)
	{
		float xP,xM,yP,yM,zP,zM;
		//Black nodes
		//laplace
		for(tz=1;tz<spltGSWBC_Z-1;tz++)
			for(ty=1;ty<spltGSWBC_Y-1;ty++)
				for(tx=0;tx<spltGSWBC_X-CUDAXTRAX;tx++)
		{
			i=tx+ty*pitchX+tz*pitchXY;
			
			float xP0,xP3,xP5,xP6;
			float yP0,yP3,yP5,yP6;
			float zP0,zP3,zP5,zP6;
	
			float shP;
	
	
	//do over P1
			shP=d_P1[i];
			
	
			zP5=shP+d_P1[i+pitchXY];
			yP3=shP+d_P1[i+pitchX];
			xP0=shP+d_P1[i-1];
	
	//do over P2
			shP=d_P2[i];
			
	
			zP6=shP+d_P2[i+pitchXY];
			yP0=shP+d_P2[i-pitchX];
			xP3=shP+d_P2[i+1];
	
	//do over P4
			shP=d_P4[i];
			
			zP0=shP+d_P4[i-pitchXY];
			yP6=shP+d_P4[i+pitchX];
			xP5=shP+d_P4[i+1];
	
	//do over P7
			shP=d_P7[i];
			
			zP3=shP+d_P7[i-pitchXY];
			yP5=shP+d_P7[i-pitchX];
			xP6=shP+d_P7[i-1];
	
	
			d_P0[i]=om1*d_P0[i]+om2d6*((xP0+yP0)+zP0);
			d_P3[i]=om1*d_P3[i]+om2d6*((xP3+yP3)+zP3);
			d_P5[i]=om1*d_P5[i]+om2d6*((xP5+yP5)+zP5);
			d_P6[i]=om1*d_P6[i]+om2d6*((xP6+yP6)+zP6);
		}
		//charges
		for(t=0;t<d_Qnum[0];t++)
		{
			i=d_Qpos[0][t];
			d_P0[i]=d_P0[i]+d_Q[0][t];
		}
		for(t=0;t<d_Qnum[3];t++)
		{
			i=d_Qpos[3][t];
			d_P3[i]=d_P3[i]+d_Q[3][t];
		}
		for(t=0;t<d_Qnum[5];t++)
		{
			i=d_Qpos[5][t];
			d_P5[i]=d_P5[i]+d_Q[5][t];
		}
		for(t=0;t<d_Qnum[6];t++)
		{
			i=d_Qpos[6][t];
			d_P6[i]=d_P6[i]+d_Q[6][t];
		}
		//Diel Border
		for(t=0;t<dc_DBnum[0];t++)
		{
			i=dc_DBpos[0][t];
		//P0 x
			xP=FMUL(dc_DielMult[PlusX][t],d_P1[i]);
			xM=FMUL(dc_DielMult[MinusX][t],d_P1[i-1]);
		//P0 y 
			yP=FMUL(dc_DielMult[PlusY][t],d_P2[i]);
			yM=FMUL(dc_DielMult[MinusY][t],d_P2[i-pitchX]);
		//P0 z
			zP=FMUL(dc_DielMult[PlusZ][t],d_P4[i]);
			zM=FMUL(dc_DielMult[MinusZ][t],d_P4[i-pitchXY]);
		
			xP=FADD(xM,xP);
			yP=FADD(yM,yP);
			zP=FADD(zM,zP);
			d_P0[i]=FADD(FADD(xP,yP),FADD(zP,d_P0[i]));
		}
		for(t=0;t<dc_DBnum[3];t++)
		{
			i=dc_DBpos[3][t];
		//P3 x
			xM=FMUL(dc_DielMult[3*6+MinusX][t],d_P2[i]);
			xP=FMUL(dc_DielMult[3*6+PlusX][t],d_P2[i+1]);
		//P3 y
			yM=FMUL(dc_DielMult[3*6+MinusY][t],d_P1[i]);
			yP=FMUL(dc_DielMult[3*6+PlusY][t],d_P1[i+pitchX]);
		//P3 z
			zP=FMUL(dc_DielMult[3*6+PlusZ][t],d_P7[i]);
			zM=FMUL(dc_DielMult[3*6+MinusZ][t],d_P7[i-pitchXY]);
		
			xP=FADD(xM,xP);
			yP=FADD(yM,yP);
			zP=FADD(zM,zP);
			d_P3[i]=FADD(FADD(xP,yP),FADD(zP,d_P3[i]));
		}
		for(t=0;t<dc_DBnum[5];t++)
		{
			i=dc_DBpos[5][t];
		//P5 x
			xM=FMUL(dc_DielMult[5*6+MinusX][t],d_P4[i]);
			xP=FMUL(dc_DielMult[5*6+PlusX][t],d_P4[i+1]);
		//P5 y
			yP=FMUL(dc_DielMult[5*6+PlusY][t],d_P7[i]);
			yM=FMUL(dc_DielMult[5*6+MinusY][t],d_P7[i-pitchX]);
		//P5 z
			zM=FMUL(dc_DielMult[5*6+MinusZ][t],d_P1[i]);
			zP=FMUL(dc_DielMult[5*6+PlusZ][t],d_P1[i+pitchXY]);
		
			xP=FADD(xM,xP);
			yP=FADD(yM,yP);
			zP=FADD(zM,zP);
			d_P5[i]=FADD(FADD(xP,yP),FADD(zP,d_P5[i]));
		}
		for(t=0;t<dc_DBnum[6];t++)
		{
			i=dc_DBpos[6][t];
			//P6 x
			xP=FMUL(dc_DielMult[6*6+PlusX][t],d_P7[i]);
			xM=FMUL(dc_DielMult[6*6+MinusX][t],d_P7[i-1]);
		//P6 y
			yM=FMUL(dc_DielMult[6*6+MinusY][t],d_P4[i]);
			yP=FMUL(dc_DielMult[6*6+PlusY][t],d_P4[i+pitchX]);
		//P6 z
			zM=FMUL(dc_DielMult[6*6+MinusZ][t],d_P2[i]);
			zP=FMUL(dc_DielMult[6*6+PlusZ][t],d_P2[i+pitchXY]);
		
			xP=FADD(xM,xP);
			yP=FADD(yM,yP);
			zP=FADD(zM,zP);
			d_P6[i]=FADD(FADD(xP,yP),FADD(zP,d_P6[i]));
		}
		//white nodes
		//laplace
		for(tz=1;tz<spltGSWBC_Z-1;tz++)
			for(ty=1;ty<spltGSWBC_Y-1;ty++)
				for(tx=0;tx<spltGSWBC_X-CUDAXTRAX;tx++)
		{
			i=tx+ty*pitchX+tz*pitchXY;
			
			float xP1,xP2,xP4,xP7;
			float yP1,yP2,yP4,yP7;
			float zP1,zP2,zP4,zP7;
	
			float shP;
	
	//do dc_P[0]
			shP=d_P0[i];
			
			zP4=shP+d_P0[i+pitchXY];
			yP2=shP+d_P0[i+pitchX];
			xP1=shP+d_P0[i+1];
	
	//do d_P[3]
			shP=d_P3[i];
			
			zP7=shP+d_P3[i+pitchXY];
			yP1=shP+d_P3[i-pitchX];
			xP2=shP+d_P3[i-1];
	
	//do d_P[5]
			shP=d_P5[i];
			
			zP1=shP+d_P5[i-pitchXY];
			yP7=shP+d_P5[i+pitchX];
			xP4=shP+d_P5[i-1];
	
	//do d_P6
			shP=d_P6[i];
			
			zP2=shP+d_P6[i-pitchXY];
			yP4=shP+d_P6[i-pitchX];
			xP7=shP+d_P6[i+1];
	
	//d_P6[i]=FADD(FMUL(om1,d_P6[i]),FMUL(om2d6,FADD(FADD(xP6,yP6),zP6)));
	
			d_P1[i]=om1*d_P1[i]+om2d6*((xP1+yP1)+zP1);
			d_P2[i]=om1*d_P2[i]+om2d6*((xP2+yP2)+zP2);
			d_P4[i]=om1*d_P4[i]+om2d6*((xP4+yP4)+zP4);
			d_P7[i]=om1*d_P7[i]+om2d6*((xP7+yP7)+zP7);
		}
		//charges
		for(t=0;t<d_Qnum[1];t++)
		{
			i=d_Qpos[1][t];
			d_P1[i]=d_P1[i]+d_Q[1][t];
		}
		for(t=0;t<d_Qnum[2];t++)
		{
			i=d_Qpos[2][t];
			d_P2[i]=d_P2[i]+d_Q[2][t];
		}
		for(t=0;t<d_Qnum[4];t++)
		{
			i=d_Qpos[4][t];
			d_P4[i]=d_P4[i]+d_Q[4][t];
		}
		for(t=0;t<d_Qnum[7];t++)
		{
			i=d_Qpos[7][t];
			d_P7[i]=d_P7[i]+d_Q[7][t];
		}
		//Diel Border
		for(t=0;t<dc_DBnum[1];t++)
		{
			i=dc_DBpos[1][t];
		//P1 x
			xM=FMUL(dc_DielMult[1*6+MinusX][t],d_P0[i]);
			xP=FMUL(dc_DielMult[1*6+PlusX][t],d_P0[i+1]);
		//P1 y
			yP=FMUL(dc_DielMult[1*6+PlusY][t],d_P3[i]);
			yM=FMUL(dc_DielMult[1*6+MinusY][t],d_P3[i-pitchX]);
		//P1 z
			zP=FMUL(dc_DielMult[1*6+PlusZ][t],d_P5[i]);
			zM=FMUL(dc_DielMult[1*6+MinusZ][t],d_P5[i-pitchXY]);
		
			xP=FADD(xM,xP);
			yP=FADD(yM,yP);
			zP=FADD(zM,zP);
			d_P1[i]=FADD(FADD(xP,yP),FADD(zP,d_P1[i]));
		}
		for(t=0;t<dc_DBnum[2];t++)
		{
			i=dc_DBpos[2][t];
		//P2 x
			xP=FMUL(dc_DielMult[2*6+PlusX][t],d_P3[i]);
			xM=FMUL(dc_DielMult[2*6+MinusX][t],d_P3[i-1]);
		//P2 y
			yM=FMUL(dc_DielMult[2*6+MinusY][t],d_P0[i]);
			yP=FMUL(dc_DielMult[2*6+PlusY][t],d_P0[i+pitchX]);
		//P2 z
			zP=FMUL(dc_DielMult[2*6+PlusZ][t],d_P6[i]);
			zM=FMUL(dc_DielMult[2*6+MinusZ][t],d_P6[i-pitchXY]);
		
			xP=FADD(xM,xP);
			yP=FADD(yM,yP);
			zP=FADD(zM,zP);
			d_P2[i]=FADD(FADD(xP,yP),FADD(zP,d_P2[i]));
		}
		for(t=0;t<dc_DBnum[4];t++)
		{
			i=dc_DBpos[4][t];
		//P4x
			xP=FMUL(dc_DielMult[4*6+PlusX][t],d_P5[i]);
			xM=FMUL(dc_DielMult[4*6+MinusX][t],d_P5[i-1]);
		//P4 y
			yP=FMUL(dc_DielMult[4*6+PlusY][t],d_P6[i]);
			yM=FMUL(dc_DielMult[4*6+MinusY][t],d_P6[i-pitchX]);
		//P4 z
			zM=FMUL(dc_DielMult[4*6+MinusZ][t],d_P0[i]);
			zP=FMUL(dc_DielMult[4*6+PlusZ][t],d_P0[i+pitchXY]);
		
			xP=FADD(xM,xP);
			yP=FADD(yM,yP);
			zP=FADD(zM,zP);
			d_P4[i]=FADD(FADD(xP,yP),FADD(zP,d_P4[i]));
		}
		for(t=0;t<dc_DBnum[7];t++)
		{
			i=dc_DBpos[7][t];
		//P7 x
			xM=FMUL(dc_DielMult[7*6+MinusX][t],d_P6[i]);
			xP=FMUL(dc_DielMult[7*6+PlusX][t],d_P6[i+1]);
		//P7 y
			yM=FMUL(dc_DielMult[7*6+MinusY][t],d_P5[i]);
			yP=FMUL(dc_DielMult[7*6+PlusY][t],d_P5[i+pitchX]);
		//P7 z
			zM=FMUL(dc_DielMult[7*6+MinusZ][t],d_P3[i]);
			zP=FMUL(dc_DielMult[7*6+PlusZ][t],d_P3[i+pitchXY]);
		
			xP=FADD(xM,xP);
			yP=FADD(yM,yP);
			zP=FADD(zM,zP);
			d_P7[i]=FADD(FADD(xP,yP),FADD(zP,d_P7[i]));
		}
		
		if(iteration%ConvergenceCheck==0)
		{
			double OldTotalEnergy=totalEnergy;
			totalEnergy=0.0;
			EnergyCharge=0.0;
			for(i=0;i<8;i++)
				for(j=0;j<PS.Qnum[i];j++)
			{
					dtmp1=double(PS.P[i][PS.Qpos[i][j]])*double(PS.Q[i][j])/double(PS.Qmult[i][j]);
					//ftmp=PS.P[i][PS.Qpos[i][j]]*PS.Q[i][j]/PS.Qmult[i][j];
					EnergyCharge+=dtmp1;
			}
			totalEnergy=EnergyCharge/(fpoh*2.0);
			
			totalChange=totalEnergy-OldTotalEnergy;
			relativeChange=totalChange/totalEnergy;
			
			pnpPrintGroup0("<PoissonIterations Nit=\"%8d\" E=\"%20.16e\" dE=\"%.4e\" rel.E=\"%.4e\" ConvFac=\"%.4e\"/>\n", iteration, totalEnergy, totalChange, relativeChange,ConvFac);
			if(PS.Tolerance!=0.0)
			{
				if(fabs(relativeChange)<=PS.Tolerance)
				{
					printf("Solver has reached the requiered tolerance level\n");
					break;
				}
			}
		}
	}
	m_ContWorld->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////
PoissonSolverOnCudaDouble::PoissonSolverOnCudaDouble()
	: GenericSolver()
{
	InitZero();
}


PoissonSolverOnCudaDouble::~PoissonSolverOnCudaDouble()
{
	Clear();
}
int PoissonSolverOnCudaDouble::InitZero()
{
	HaObject::SetName("PoissonSolverOnCudaDouble");
	SolverStr.push_back("Auto");
	SolverStr.push_back("CPU");
	SolverStr.push_back("GPU");
	
	CudaParm.BS_X=96;
	CudaParm.BS_Y=2;
	CudaParm.BS_Z=1;
	CudaParm.BS_XY=CudaParm.BS_X*CudaParm.BS_Y;
	CudaParm.BS_XYZ=CudaParm.BS_XY*CudaParm.BS_Z;
	
	CudaParm.Qblock=32;
	CudaParm.DBblock=512;
	m_ContWorld=NULL;
	
	int i,j;
	PS.Relaxation=1.6;
	PS.MaxIterations=2000;
	PS.ConvergenceCheck=20;
	PS.Tolerance=0.0;
	PS.GS[0]=1;
	PS.GS[1]=1;
	PS.GS[2]=1;
	PS.spltGSWBC[0]=1;
	PS.spltGSWBC[1]=1;
	PS.spltGSWBC[2]=1;
	PS.GridScale=1.0;
	PS.spltGSWBC_X=1;
	PS.spltGSWBC_XY=1;
	PS.spltGSWBC_XYZ=1;
	
	PS.TotalEnergy=0.0;
	PS.AvrOverChecks=0;
	PS.TEavr=0.0;
	PS.stdevTE=0.0;
	for(i=0;i<8;i++)
	{
		PS.P[i]=NULL;
		PS.Qnum[i]=0;
		PS.Q[i]=NULL;
		PS.Qmult[i]=NULL;
		PS.Qpos[i]=NULL;
		
		PS.DielBordNum[i]=0;
		PS.DielBordPos[i]=NULL;
		for(j=0;j<6;j++)
			PS.DielMult[i][j]=NULL;
	}
	return EXIT_SUCCESS;
}
int PoissonSolverOnCudaDouble::Clear()
{
	int i,j;
	for(i=0;i<8;i++)
	{
		DeleteCArray(PS.P[i]);
		DeleteCArray(PS.Q[i]);
		DeleteCArray(PS.Qmult[i]);
		DeleteCArray(PS.Qpos[i]);
		if(PS.DielBordPos[i]!=NULL)
		{
			delete [] PS.DielBordPos[i];
			PS.DielBordPos[i]=NULL;
		}
		//DeleteCArray(PS.DielBordPos[i]);
		for(j=0;j<6;j++)
			DeleteCArray(PS.DielMult[i][j]);
	}
	return EXIT_SUCCESS;
}
int PoissonSolverOnCudaDouble::SaveXML(TiXmlElement* Elt, HaContext* p_ctxt )
{
	return EXIT_SUCCESS;
}
int PoissonSolverOnCudaDouble::LoadXML(const TiXmlElement* Elt, HaContext* p_ctxt )
{
	if(strncmp(HaObject::GetCStrName(),Elt->Value(),13))
	{
		fprintf(stderr,"ERROR: Wrong XML Element %s, expecting %s\n",Elt->Value(),HaObject::GetCStrName());
		return EXIT_FAILURE;
	}
	//Read Primary Parameters
	Elt->GetIntAttribute("MaxIterations",&PS.MaxIterations);
	Elt->GetDoubleAttribute("Relaxation",&PS.Relaxation);
	Elt->GetIntAttribute("ConvergenceCheck",&PS.ConvergenceCheck);
	
	if(Elt->GetStdStrIndex("Solver",&solver,SolverStr)!=EXIT_SUCCESS)solver=2;
	if(solver==0)solver=2;
	return EXIT_SUCCESS;
}
int PoissonSolverOnCudaDouble::SetContWorld(ContWorld* _world)
{
	m_ContWorld=_world;
	PS.GS[0]=m_ContWorld->GridSizeGlobal[0]-2;
	PS.GS[1]=m_ContWorld->GridSizeGlobal[1]-2;
	PS.GS[2]=m_ContWorld->GridSizeGlobal[2]-2;
	PS.spltGSWBC[0]=PS.GS[0]/2+CUDAXTRAX;
	PS.spltGSWBC[1]=PS.GS[1]/2+2;
	PS.spltGSWBC[2]=PS.GS[2]/2+2;
	PS.spltGSWBC_X=PS.spltGSWBC[0];
	PS.spltGSWBC_XY=PS.spltGSWBC[0]*PS.spltGSWBC[1];
	PS.spltGSWBC_XYZ=PS.spltGSWBC[0]*PS.spltGSWBC[1]*PS.spltGSWBC[2];
	PS.GridScale=m_ContWorld->GridScale;
	return EXIT_SUCCESS;
}
int PoissonSolverOnCudaDouble::ShowParameters()
{
	return EXIT_SUCCESS;
}
int PoissonSolverOnCudaDouble::ShowProperties()
{
	return EXIT_SUCCESS;
}
int PoissonSolverOnCudaDouble::InitSolver()
{
	DbgPrint0("PoissonSolverOnCudaDouble::InitSolver()\n");
	int i,j,k,GrdPnt;
	int ix,iy,iz;
	int t,tx,ty,tz;
	double ftmp1,ftmp2,ftmp3;
	double om2=PS.Relaxation;
	double om2d6=PS.Relaxation/6.0;
	NodeIndexing* NIndexing=m_ContWorld->NIndexing;
	NodeIndex* NIndex=NIndexing->NIndex;
	double Eps[NodeIndexMaxValues];
	for(i=0;i<NodeIndexMaxValues;i++)
		Eps[i]=double(NIndexing->Eps[i]);
	
	
	NodeIndex specChargeMask=NodeIndexing::ChargeMask;
	NodeIndex DielBoarderMask=NodeIndexing::DielBoarderMask;
	
	int spltGSWBC_XYZ=PS.spltGSWBC[0]*PS.spltGSWBC[1]*PS.spltGSWBC[2];
	//allocate host split potential
	for(i=0;i<8;i++)
	{
		PS.P[i]=new double[spltGSWBC_XYZ];
		for(j=0;j<spltGSWBC_XYZ;j++)
			PS.P[i][j]=0.0f;
	}
	
	
	int pGS_X=m_ContWorld->GS_X;
	int pGS_Y=m_ContWorld->GS_Y;
	int pGS_Z=m_ContWorld->GS_Z;
	int pGS_XY=pGS_X*pGS_Y;
	int pGS_XYZ=pGS_X*pGS_Y*pGS_Z;
	
	int izgrid,iygrid;
	//cound Q and BielBorder
	for(i=0;i<8;i++)
	{
		PS.Qnum[i]=0;
		PS.DielBordNum[i]=0;
	}
	int QCount=0;
	int DBCount=0;
	for(iz=1;iz<pGS_Z-1;iz++)
	{
		izgrid = iz*pGS_XY;
		for(iy=1;iy<pGS_Y-1;iy++)
		{
			iygrid = izgrid+iy*pGS_X;
			for(ix=1;ix<pGS_X-1;ix++)
			{
				GrdPnt = iygrid+ix;
				if(NIndex[GrdPnt]&specChargeMask)
				{
					i=((ix - 1)%2)+((iy + 1)%2)*2+((iz + 1)%2)*4;
					PS.Qnum[i]++;
					QCount++;
				}
				if(NIndex[GrdPnt]&DielBoarderMask)
				{
					i=((ix - 1)%2)+((iy + 1)%2)*2+((iz + 1)%2)*4;
					PS.DielBordNum[i]++;
					DBCount++;
				}
			}
		}
	}
	for(i=0;i<8;i++)
		DbgPrint0("PS.Qnum[%d]=%d\n",i,PS.Qnum[i]);
	DbgPrint0("QCount %d\n",QCount);
	for(i=0;i<8;i++)
		DbgPrint0("PS.DielBordNum[%d]=%d\n",i,PS.DielBordNum[i]);
	DbgPrint0("DBCount %d\n",DBCount);
	//allocate split charge and diel border
	int iQcount[8];
	int iDBcount[8];
	for(i=0;i<8;i++)
	{
		if(PS.Qnum[i]>0)
		{
			PS.Q[i]=new double[PS.Qnum[i]];
			PS.Qpos[i]=new int[PS.Qnum[i]];
			PS.Qmult[i]=new double[PS.Qnum[i]];
			for(j=0;j<PS.Qnum[i];j++)
			{
				PS.Q[i][j]=0.0f;
				PS.Qpos[i][j]=-1;
				PS.Qmult[i][j]=0.0f;
			}
		}
		else
		{
			PS.Q[i]=NULL;
			PS.Qpos[i]=NULL;
			PS.Qmult[i]=NULL;
		}
		if(PS.DielBordNum[i]>0)
		{
			PS.DielBordPos[i]=new int[PS.DielBordNum[i]];
			for(j=0;j<PS.DielBordNum[i];j++)
				PS.DielBordPos[i][j]=-1;
			for(k=0;k<6;k++)
			{
				PS.DielMult[i][k]=new double[PS.DielBordNum[i]];
				for(j=0;j<PS.DielBordNum[i];j++)
					PS.DielMult[i][k][j]=0.0;
			}
		}
		else
		{
			PS.DielBordPos[i]==NULL;
			for(k=0;k<6;k++)
				PS.DielMult[i][k]=NULL;
		}
		iQcount[i]=0;
		iDBcount[i]=0;
	}
	//fill split charge and diel border
	QCount=0;
	for(iz=1;iz<pGS_Z-1;iz++)
	{
		izgrid = iz*pGS_XY;
		for(iy=1;iy<pGS_Y-1;iy++)
		{
			iygrid = izgrid+iy*pGS_X;
			for(ix=1;ix<pGS_X-1;ix++)
			{
				GrdPnt = iygrid+ix;
				if(NIndex[GrdPnt]&specChargeMask)
				{
					tx=(ix - 1)/2;
					ty=(iy + 1)/2;
					tz=(iz + 1)/2;
					t=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
					i=((ix - 1)%2)+((iy + 1)%2)*2+((iz + 1)%2)*4;
					
					PS.Qpos[i][iQcount[i]]=t;
					if(NIndex[GrdPnt]&DielBoarderMask)
					{
						ftmp1=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft] + Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft] + Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft] + Eps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft] + Eps[(NIndex[GrdPnt-pGS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft] + Eps[(NIndex[GrdPnt-pGS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
						PS.Qmult[i][iQcount[i]]=om2/ftmp1;
						PS.Q[i][iQcount[i]]=PS.Qmult[i][iQcount[i]]*NIndexing->Q[QCount];
					}
					else
					{
						PS.Qmult[i][iQcount[i]]=om2d6/Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
						PS.Q[i][iQcount[i]]=PS.Qmult[i][iQcount[i]]*NIndexing->Q[QCount];
					}
					
					
					iQcount[i]++;
					QCount++;
				}
				if(NIndex[GrdPnt]&DielBoarderMask)
				{
					tx=(ix - 1)/2;
					ty=(iy + 1)/2;
					tz=(iz + 1)/2;
					t=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
					i=((ix - 1)%2)+((iy + 1)%2)*2+((iz + 1)%2)*4;
					
					PS.DielBordPos[i][iDBcount[i]]=t;
					
					PS.DielMult[i][PlusX][iDBcount[i]] =Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
					PS.DielMult[i][PlusY][iDBcount[i]] =Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
					PS.DielMult[i][PlusZ][iDBcount[i]] =Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
					PS.DielMult[i][MinusX][iDBcount[i]]=Eps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
					PS.DielMult[i][MinusY][iDBcount[i]]=Eps[(NIndex[GrdPnt-pGS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
					PS.DielMult[i][MinusZ][iDBcount[i]]=Eps[(NIndex[GrdPnt-pGS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
					
					ftmp1=PS.DielMult[i][PlusX][iDBcount[i]]+PS.DielMult[i][MinusX][iDBcount[i]]+PS.DielMult[i][PlusY][iDBcount[i]]+PS.DielMult[i][MinusY][iDBcount[i]]+PS.DielMult[i][PlusZ][iDBcount[i]]+PS.DielMult[i][MinusZ][iDBcount[i]];
					ftmp1=PS.Relaxation/ftmp1;
					
					for(k=0;k<6;k++)
					{
						PS.DielMult[i][k][iDBcount[i]] = ftmp1*PS.DielMult[i][k][iDBcount[i]]-om2d6;
					}
					iDBcount[i]++;
				}
			}
		}
	}
	return EXIT_SUCCESS;
}
extern "C" int DoPSolverOnCudaDouble(PoissonSolverOnCudaParamStruct* CudaParm,PSolverOnCudaStructDouble* PS);
int PoissonSolverOnCudaDouble::Solve()
{
	DbgPrint0("PoissonSolverOnCudaDouble::Solve()\n");
	int i;
	if(m_ContWorld->Potential == NULL)
	{
		m_ContWorld->Potential = new float[m_ContWorld->GS_XYZ];
		for(i=0;i<m_ContWorld->GS_XYZ;i++)
			m_ContWorld->Potential[i]= 0.0;
	}
	
	
	CopyPotBC2SplitPot();
	CopyIntPot2SplitPot();printf("dimBlock [%d,%d,%d]\n",CudaParm.BS_X,CudaParm.BS_Y,CudaParm.BS_Z);
	if(solver==1)
		DoPoissonSolverCudaOnCPU();
	else
	{
		int halfGS_X=(m_ContWorld->GS_X-2)/2;
		int BS_X;
		for(CudaParm.BS_X=96;CudaParm.BS_X>=32;CudaParm.BS_X=CudaParm.BS_X-32)
		{
			DbgPrint0("halfGS_XCudaParm.BS_X %d %d %d %d\n",m_ContWorld->GS_X,halfGS_X,CudaParm.BS_X,halfGS_X%CudaParm.BS_X);
			if(halfGS_X%CudaParm.BS_X==0)
			{
				DbgPrint0("BS_X=%d and it is optimal\n",CudaParm.BS_X);
				break;
			}
		}
		if(halfGS_X%32!=0)
		{
			if(halfGS_X<96)
			{
				CudaParm.BS_X=halfGS_X;
				DbgPrint0("BS_X=%d and it is NOT optimal\n",CudaParm.BS_X);
			}
			else
			{
				pnpError("Can not find block size for GPU, change GridSize of the system\n");
			}
		}
		CudaParm.BS_Y=2;
		CudaParm.BS_Z=1;
		CudaParm.BS_XY=CudaParm.BS_X*CudaParm.BS_Y;
		CudaParm.BS_XYZ=CudaParm.BS_XY*CudaParm.BS_Z;
#if defined(WITH_CUDA)
		DoPSolverOnCudaDouble(&CudaParm,&PS);
#endif
	}
	CopySplitPot2IntPot();
	
	
// 	CopyPotBC2SplitPot();
// 	int Qblocks[6]={16,32,64,128,256,512};
// 	for(int i=0;i<6;i++)
// 	{
// 		CudaParm.Qblock=Qblocks[i];
// 		CopyIntPot2SplitPot();
// 		DoPoissonSolverOnCuda(PS,CudaParm);
// 	}
	//calculate energy
	CalcSystemEnergy();
	pnpPrintGroup0("<PoissonFinal E=\"%.10e\"/>\n", m_ContWorld->SystemEnergy);
	return EXIT_SUCCESS;
}
int PoissonSolverOnCudaDouble::CalcSystemEnergy()
{
	int i,j;
	static double oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	double Dev=0.0;
	double tmp,SumSQ=0.0,SumAbs=0.0;
	double EnergyCharge = 0.0;
	double EnergySingular = 0.0;
	double EnergyQmob = 0.0;
	double totalEnergy=0.0;
	double fpoh = 4.0*M_PI*m_ContWorld->GridScale;
	double dphi=0.0,maxdphi=0.0;
	int GrdPntMaxDPhi=-1;
	double om2d6=PS.Relaxation/6.0;
	for(i=0;i<8;i++)
	{
		for(j=0;j<PS.Qnum[i];j++)
		{
			tmp=double(PS.P[i][PS.Qpos[i][j]])*double(PS.Q[i][j])/double(PS.Qmult[i][j]);
			EnergyCharge+=tmp;
		}
	}
	totalEnergy=EnergyCharge/(fpoh*2.0);
	m_ContWorld->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
int PoissonSolverOnCudaDouble::CopyPotBC2SplitPot()
{
	float* Pot=m_ContWorld->Potential;
	int GrdPnt0,ix,iy,iz;
	int i,tx,ty,tz;
	int PotGS_X=m_ContWorld->GS_X;
	int PotGS_XY=m_ContWorld->GS_XY;
	//copy BC Potential to split potential
	//XY
	for(ty=1;ty<PS.spltGSWBC[1]-1;ty++)
		for(tx=0;tx<PS.spltGSWBC[0]-CUDAXTRAX;tx++)
	{
		tz=0;
		ix = 2*tx+1;iy = 2*ty-1;iz = 2*tz-1;
		i=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
		GrdPnt0=ix+iy*PotGS_X+iz*PotGS_XY;
		
		PS.P[4][i]=Pot[GrdPnt0+PotGS_XY];
		PS.P[5][i]=Pot[GrdPnt0+PotGS_XY+1];
		PS.P[6][i]=Pot[GrdPnt0+PotGS_XY+PotGS_X];
		PS.P[7][i]=Pot[GrdPnt0+PotGS_XY+1+PotGS_X];
			
		tz=PS.spltGSWBC[2]-1;
		ix = 2*tx+1;iy = 2*ty-1;iz = 2*tz-1;	
		i=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
		GrdPnt0=ix+iy*PotGS_X+iz*PotGS_XY;
		
		PS.P[0][i]=Pot[GrdPnt0];
		PS.P[1][i]=Pot[GrdPnt0+1];
		PS.P[2][i]=Pot[GrdPnt0+PotGS_X];
		PS.P[3][i]=Pot[GrdPnt0+1+PotGS_X];
		
	}
	//ZX
	for(tz=1;tz<PS.spltGSWBC[2]-1;tz++)
		for(tx=0;tx<PS.spltGSWBC[0]-CUDAXTRAX;tx++)
	{
		ty=0;
		ix = 2*tx+1;iy = 2*ty-1;iz = 2*tz-1;
		i=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
		GrdPnt0=ix+iy*PotGS_X+iz*PotGS_XY;
			
		PS.P[2][i]=Pot[GrdPnt0+PotGS_X];
		PS.P[3][i]=Pot[GrdPnt0+1+PotGS_X];
		PS.P[6][i]=Pot[GrdPnt0+PotGS_XY+PotGS_X];
		PS.P[7][i]=Pot[GrdPnt0+PotGS_XY+1+PotGS_X];
			
		ty=PS.spltGSWBC[1]-1;
		ix = 2*tx+1;iy = 2*ty-1;iz = 2*tz-1;
		i=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
		GrdPnt0=ix+iy*PotGS_X+iz*PotGS_XY;
			
		PS.P[0][i]=Pot[GrdPnt0];
		PS.P[1][i]=Pot[GrdPnt0+1];
		PS.P[4][i]=Pot[GrdPnt0+PotGS_XY];
		PS.P[5][i]=Pot[GrdPnt0+PotGS_XY+1];
	}
	//YZ
	for(tz=1;tz<PS.spltGSWBC[2]-1;tz++)
		for(ty=1;ty<PS.spltGSWBC[1]-1;ty++)
	{
		tx=-1;
		ix = 2*tx+1;iy = 2*ty-1;iz = 2*tz-1;
		i=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
		GrdPnt0=ix+iy*PotGS_X+iz*PotGS_XY;
		PS.P[1][i]=Pot[GrdPnt0+1];
		PS.P[3][i]=Pot[GrdPnt0+1+PotGS_X];
		PS.P[5][i]=Pot[GrdPnt0+PotGS_XY+1];
		PS.P[7][i]=Pot[GrdPnt0+PotGS_XY+1+PotGS_X];
			
		tx=PS.spltGSWBC[0]-CUDAXTRAX;
		ix = 2*tx+1;iy = 2*ty-1;iz = 2*tz-1;
		i=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
		GrdPnt0=ix+iy*PotGS_X+iz*PotGS_XY;
		PS.P[0][i]=Pot[GrdPnt0];
		PS.P[2][i]=Pot[GrdPnt0+PotGS_X];
		PS.P[4][i]=Pot[GrdPnt0+PotGS_XY];
		PS.P[6][i]=Pot[GrdPnt0+PotGS_XY+PotGS_X];
			
	}
	return EXIT_SUCCESS;
}
int PoissonSolverOnCudaDouble::CopyIntPot2SplitPot()
{
	float* Pot=m_ContWorld->Potential;
	int GrdPnt0,ix,iy,iz;
	int i,tx,ty,tz;
	int PotGS_X=m_ContWorld->GS_X;
	int PotGS_XY=m_ContWorld->GS_XY;
	//copy internal Potential to split potential
	for(tz=1;tz<PS.spltGSWBC[2]-1;tz++)
		for(ty=1;ty<PS.spltGSWBC[1]-1;ty++)
			for(tx=0;tx<PS.spltGSWBC[0]-CUDAXTRAX;tx++)
	{
		ix = 2*tx+1;
		iy = 2*ty-1;
		iz = 2*tz-1;
			
		i=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
		GrdPnt0=ix+iy*PotGS_X+iz*PotGS_XY;
		PS.P[0][i]=Pot[GrdPnt0];
		PS.P[1][i]=Pot[GrdPnt0+1];
		PS.P[2][i]=Pot[GrdPnt0+PotGS_X];
		PS.P[3][i]=Pot[GrdPnt0+1+PotGS_X];
		PS.P[4][i]=Pot[GrdPnt0+PotGS_XY];
		PS.P[5][i]=Pot[GrdPnt0+PotGS_XY+1];
		PS.P[6][i]=Pot[GrdPnt0+PotGS_XY+PotGS_X];
		PS.P[7][i]=Pot[GrdPnt0+PotGS_XY+1+PotGS_X];
	}
	return EXIT_SUCCESS;
}
int PoissonSolverOnCudaDouble::CopySplitPot2IntPot()
{
	float* Pot=m_ContWorld->Potential;
	int GrdPnt0,ix,iy,iz;
	int i,tx,ty,tz;
	int PotGS_X=m_ContWorld->GS_X;
	int PotGS_XY=m_ContWorld->GS_XY;
	//copy internal part of split potential to internal part of Potential
	for(tz=1;tz<PS.spltGSWBC[2]-1;tz++)
		for(ty=1;ty<PS.spltGSWBC[1]-1;ty++)
			for(tx=0;tx<PS.spltGSWBC[0]-CUDAXTRAX;tx++)
	{
		ix = 2*tx+1;
		iy = 2*ty-1;
		iz = 2*tz-1;
			
		i=tx+ty*PS.spltGSWBC_X+tz*PS.spltGSWBC_XY;
		GrdPnt0=ix+iy*PotGS_X+iz*PotGS_XY;
		Pot[GrdPnt0]=PS.P[0][i];
		Pot[GrdPnt0+1]=PS.P[1][i];
		Pot[GrdPnt0+PotGS_X]=PS.P[2][i];
		Pot[GrdPnt0+1+PotGS_X]=PS.P[3][i];
		Pot[GrdPnt0+PotGS_XY]=PS.P[4][i];
		Pot[GrdPnt0+PotGS_XY+1]=PS.P[5][i];
		Pot[GrdPnt0+PotGS_XY+PotGS_X]=PS.P[6][i];
		Pot[GrdPnt0+PotGS_XY+1+PotGS_X]=PS.P[7][i];
	}
	return EXIT_SUCCESS;
}

int PoissonSolverOnCudaDouble::DoPoissonSolverCudaOnCPU()
{
	pnpPrintGroup0("Solving Poisson Equation using algorithm made for CUDA\n");
	int i,j,k,t,tx,ty,tz;
	int pitchX=PS.spltGSWBC_X;
	int pitchXY=PS.spltGSWBC_XY;
	int spltGSWBC_X=PS.spltGSWBC[0];
	int spltGSWBC_Y=PS.spltGSWBC[1];
	int spltGSWBC_Z=PS.spltGSWBC[2];
	int MaxIterations=PS.MaxIterations;
	int ConvergenceCheck=PS.ConvergenceCheck;
	
	
	int iteration;
	double om2 = PS.Relaxation;
	double om1 = 1.0-om2;
	double om2d6 = om2/6.0;
	
	double *d_P0=PS.P[0];
	double *d_P1=PS.P[1];
	double *d_P2=PS.P[2];
	double *d_P3=PS.P[3];
	double *d_P4=PS.P[4];
	double *d_P5=PS.P[5];
	double *d_P6=PS.P[6];
	double *d_P7=PS.P[7];
	
	double ftmp1,ftmp2;
	double dtmp1,dtmp2;
	
	double fpoh = 4.0*M_PI*m_ContWorld->GridScale;
	double ftmp;
	double totalEnergy=0.0;
	double totalChange=0.0;
	double relativeChange=0.0;
	double ConvFac=0.0;
	double EnergyCharge = 0.0;
	
	int d_Qnum[8];
	double *d_Q[8];
	int *d_Qpos[8];
	
	double* dc_DielMult[48];
	int dc_DBnum[8];
	int* dc_DBpos[8];
	for(i=0;i<8;i++)
	{
		dc_DBpos[i]=NULL;
		for(k=0;k<6;k++)
			dc_DielMult[i*6+k]=NULL;
		if(PS.DielBordNum[i]>0)
		{
			for(k=0;k<6;k++)
			{
				dc_DielMult[i*6+k]=PS.DielMult[i][k];
			}
			dc_DBpos[i]=PS.DielBordPos[i];
		}
		dc_DBnum[i]=PS.DielBordNum[i];
	}
	
	
	
	for(i=0;i<8;i++)
	{
		d_Qnum[i]=PS.Qnum[i];
		d_Q[i]=PS.Q[i];
		d_Qpos[i]=PS.Qpos[i];
	}
	for(iteration=1;iteration<=MaxIterations;iteration++)
	{
		float xP,xM,yP,yM,zP,zM;
		//Black nodes
		//laplace
		for(tz=1;tz<spltGSWBC_Z-1;tz++)
			for(ty=1;ty<spltGSWBC_Y-1;ty++)
				for(tx=0;tx<spltGSWBC_X-CUDAXTRAX;tx++)
		{
			i=tx+ty*pitchX+tz*pitchXY;
			
			double xP0,xP3,xP5,xP6;
			double yP0,yP3,yP5,yP6;
			double zP0,zP3,zP5,zP6;
	
			double shP;
	
	
	//do over P1
			shP=d_P1[i];
			
	
			zP5=shP+d_P1[i+pitchXY];
			yP3=shP+d_P1[i+pitchX];
			xP0=shP+d_P1[i-1];
	
	//do over P2
			shP=d_P2[i];
			
	
			zP6=shP+d_P2[i+pitchXY];
			yP0=shP+d_P2[i-pitchX];
			xP3=shP+d_P2[i+1];
	
	//do over P4
			shP=d_P4[i];
			
			zP0=shP+d_P4[i-pitchXY];
			yP6=shP+d_P4[i+pitchX];
			xP5=shP+d_P4[i+1];
	
	//do over P7
			shP=d_P7[i];
			
			zP3=shP+d_P7[i-pitchXY];
			yP5=shP+d_P7[i-pitchX];
			xP6=shP+d_P7[i-1];
	
	
			d_P0[i]=om1*d_P0[i]+om2d6*((xP0+yP0)+zP0);
			d_P3[i]=om1*d_P3[i]+om2d6*((xP3+yP3)+zP3);
			d_P5[i]=om1*d_P5[i]+om2d6*((xP5+yP5)+zP5);
			d_P6[i]=om1*d_P6[i]+om2d6*((xP6+yP6)+zP6);
		}
		//charges
		for(t=0;t<d_Qnum[0];t++)
		{
			i=d_Qpos[0][t];
			PS.P[0][i]=PS.P[0][i]+d_Q[0][t];
		}
		for(t=0;t<d_Qnum[3];t++)
		{
			i=d_Qpos[3][t];
			PS.P[3][i]=PS.P[3][i]+d_Q[3][t];
		}
		for(t=0;t<d_Qnum[5];t++)
		{
			i=d_Qpos[5][t];
			PS.P[5][i]=PS.P[5][i]+d_Q[5][t];
		}
		for(t=0;t<d_Qnum[6];t++)
		{
			i=d_Qpos[6][t];
			PS.P[6][i]=PS.P[6][i]+d_Q[6][t];
		}
		//Diel Border
		for(t=0;t<dc_DBnum[0];t++)
		{
			i=dc_DBpos[0][t];
		//P0 x
			xP=FMUL(dc_DielMult[PlusX][t],d_P1[i]);
			xM=FMUL(dc_DielMult[MinusX][t],d_P1[i-1]);
		//P0 y 
			yP=FMUL(dc_DielMult[PlusY][t],d_P2[i]);
			yM=FMUL(dc_DielMult[MinusY][t],d_P2[i-pitchX]);
		//P0 z
			zP=FMUL(dc_DielMult[PlusZ][t],d_P4[i]);
			zM=FMUL(dc_DielMult[MinusZ][t],d_P4[i-pitchXY]);
		
			xP=FADD(xM,xP);
			yP=FADD(yM,yP);
			zP=FADD(zM,zP);
			d_P0[i]=FADD(FADD(xP,yP),FADD(zP,d_P0[i]));
		}
		for(t=0;t<dc_DBnum[3];t++)
		{
			i=dc_DBpos[3][t];
		//P3 x
			xM=FMUL(dc_DielMult[3*6+MinusX][t],d_P2[i]);
			xP=FMUL(dc_DielMult[3*6+PlusX][t],d_P2[i+1]);
		//P3 y
			yM=FMUL(dc_DielMult[3*6+MinusY][t],d_P1[i]);
			yP=FMUL(dc_DielMult[3*6+PlusY][t],d_P1[i+pitchX]);
		//P3 z
			zP=FMUL(dc_DielMult[3*6+PlusZ][t],d_P7[i]);
			zM=FMUL(dc_DielMult[3*6+MinusZ][t],d_P7[i-pitchXY]);
		
			xP=FADD(xM,xP);
			yP=FADD(yM,yP);
			zP=FADD(zM,zP);
			d_P3[i]=FADD(FADD(xP,yP),FADD(zP,d_P3[i]));
		}
		for(t=0;t<dc_DBnum[5];t++)
		{
			i=dc_DBpos[5][t];
		//P5 x
			xM=FMUL(dc_DielMult[5*6+MinusX][t],d_P4[i]);
			xP=FMUL(dc_DielMult[5*6+PlusX][t],d_P4[i+1]);
		//P5 y
			yP=FMUL(dc_DielMult[5*6+PlusY][t],d_P7[i]);
			yM=FMUL(dc_DielMult[5*6+MinusY][t],d_P7[i-pitchX]);
		//P5 z
			zM=FMUL(dc_DielMult[5*6+MinusZ][t],d_P1[i]);
			zP=FMUL(dc_DielMult[5*6+PlusZ][t],d_P1[i+pitchXY]);
		
			xP=FADD(xM,xP);
			yP=FADD(yM,yP);
			zP=FADD(zM,zP);
			d_P5[i]=FADD(FADD(xP,yP),FADD(zP,d_P5[i]));
		}
		for(t=0;t<dc_DBnum[6];t++)
		{
			i=dc_DBpos[6][t];
			//P6 x
			xP=FMUL(dc_DielMult[6*6+PlusX][t],d_P7[i]);
			xM=FMUL(dc_DielMult[6*6+MinusX][t],d_P7[i-1]);
		//P6 y
			yM=FMUL(dc_DielMult[6*6+MinusY][t],d_P4[i]);
			yP=FMUL(dc_DielMult[6*6+PlusY][t],d_P4[i+pitchX]);
		//P6 z
			zM=FMUL(dc_DielMult[6*6+MinusZ][t],d_P2[i]);
			zP=FMUL(dc_DielMult[6*6+PlusZ][t],d_P2[i+pitchXY]);
		
			xP=FADD(xM,xP);
			yP=FADD(yM,yP);
			zP=FADD(zM,zP);
			d_P6[i]=FADD(FADD(xP,yP),FADD(zP,d_P6[i]));
		}
		//white nodes
		//laplace
		for(tz=1;tz<spltGSWBC_Z-1;tz++)
			for(ty=1;ty<spltGSWBC_Y-1;ty++)
				for(tx=0;tx<spltGSWBC_X-CUDAXTRAX;tx++)
		{
			i=tx+ty*pitchX+tz*pitchXY;
			
			double xP1,xP2,xP4,xP7;
			double yP1,yP2,yP4,yP7;
			double zP1,zP2,zP4,zP7;
	
			double shP;
	
	//do dc_P[0]
			shP=d_P0[i];
			
			zP4=shP+d_P0[i+pitchXY];
			yP2=shP+d_P0[i+pitchX];
			xP1=shP+d_P0[i+1];
	
	//do d_P[3]
			shP=d_P3[i];
			
			zP7=shP+d_P3[i+pitchXY];
			yP1=shP+d_P3[i-pitchX];
			xP2=shP+d_P3[i-1];
	
	//do d_P[5]
			shP=d_P5[i];
			
			zP1=shP+d_P5[i-pitchXY];
			yP7=shP+d_P5[i+pitchX];
			xP4=shP+d_P5[i-1];
	
	//do d_P6
			shP=d_P6[i];
			
			zP2=shP+d_P6[i-pitchXY];
			yP4=shP+d_P6[i-pitchX];
			xP7=shP+d_P6[i+1];
	
	//d_P6[i]=FADD(FMUL(om1,d_P6[i]),FMUL(om2d6,FADD(FADD(xP6,yP6),zP6)));
	
			d_P1[i]=om1*d_P1[i]+om2d6*((xP1+yP1)+zP1);
			d_P2[i]=om1*d_P2[i]+om2d6*((xP2+yP2)+zP2);
			d_P4[i]=om1*d_P4[i]+om2d6*((xP4+yP4)+zP4);
			d_P7[i]=om1*d_P7[i]+om2d6*((xP7+yP7)+zP7);
		}
		//charges
		for(t=0;t<d_Qnum[1];t++)
		{
			i=d_Qpos[1][t];
			PS.P[1][i]=PS.P[1][i]+d_Q[1][t];
		}
		for(t=0;t<d_Qnum[2];t++)
		{
			i=d_Qpos[2][t];
			PS.P[2][i]=PS.P[2][i]+d_Q[2][t];
		}
		for(t=0;t<d_Qnum[4];t++)
		{
			i=d_Qpos[4][t];
			PS.P[4][i]=PS.P[4][i]+d_Q[4][t];
		}
		for(t=0;t<d_Qnum[7];t++)
		{
			i=d_Qpos[7][t];
			PS.P[7][i]=PS.P[7][i]+d_Q[7][t];
		}
		//Diel Border
		for(t=0;t<dc_DBnum[1];t++)
		{
			i=dc_DBpos[1][t];
		//P1 x
			xM=FMUL(dc_DielMult[1*6+MinusX][t],d_P0[i]);
			xP=FMUL(dc_DielMult[1*6+PlusX][t],d_P0[i+1]);
		//P1 y
			yP=FMUL(dc_DielMult[1*6+PlusY][t],d_P3[i]);
			yM=FMUL(dc_DielMult[1*6+MinusY][t],d_P3[i-pitchX]);
		//P1 z
			zP=FMUL(dc_DielMult[1*6+PlusZ][t],d_P5[i]);
			zM=FMUL(dc_DielMult[1*6+MinusZ][t],d_P5[i-pitchXY]);
		
			xP=FADD(xM,xP);
			yP=FADD(yM,yP);
			zP=FADD(zM,zP);
			d_P1[i]=FADD(FADD(xP,yP),FADD(zP,d_P1[i]));
		}
		for(t=0;t<dc_DBnum[2];t++)
		{
			i=dc_DBpos[2][t];
		//P2 x
			xP=FMUL(dc_DielMult[2*6+PlusX][t],d_P3[i]);
			xM=FMUL(dc_DielMult[2*6+MinusX][t],d_P3[i-1]);
		//P2 y
			yM=FMUL(dc_DielMult[2*6+MinusY][t],d_P0[i]);
			yP=FMUL(dc_DielMult[2*6+PlusY][t],d_P0[i+pitchX]);
		//P2 z
			zP=FMUL(dc_DielMult[2*6+PlusZ][t],d_P6[i]);
			zM=FMUL(dc_DielMult[2*6+MinusZ][t],d_P6[i-pitchXY]);
		
			xP=FADD(xM,xP);
			yP=FADD(yM,yP);
			zP=FADD(zM,zP);
			d_P2[i]=FADD(FADD(xP,yP),FADD(zP,d_P2[i]));
		}
		for(t=0;t<dc_DBnum[4];t++)
		{
			i=dc_DBpos[4][t];
		//P4x
			xP=FMUL(dc_DielMult[4*6+PlusX][t],d_P5[i]);
			xM=FMUL(dc_DielMult[4*6+MinusX][t],d_P5[i-1]);
		//P4 y
			yP=FMUL(dc_DielMult[4*6+PlusY][t],d_P6[i]);
			yM=FMUL(dc_DielMult[4*6+MinusY][t],d_P6[i-pitchX]);
		//P4 z
			zM=FMUL(dc_DielMult[4*6+MinusZ][t],d_P0[i]);
			zP=FMUL(dc_DielMult[4*6+PlusZ][t],d_P0[i+pitchXY]);
		
			xP=FADD(xM,xP);
			yP=FADD(yM,yP);
			zP=FADD(zM,zP);
			d_P4[i]=FADD(FADD(xP,yP),FADD(zP,d_P4[i]));
		}
		for(t=0;t<dc_DBnum[7];t++)
		{
			i=dc_DBpos[7][t];
		//P7 x
			xM=FMUL(dc_DielMult[7*6+MinusX][t],d_P6[i]);
			xP=FMUL(dc_DielMult[7*6+PlusX][t],d_P6[i+1]);
		//P7 y
			yM=FMUL(dc_DielMult[7*6+MinusY][t],d_P5[i]);
			yP=FMUL(dc_DielMult[7*6+PlusY][t],d_P5[i+pitchX]);
		//P7 z
			zM=FMUL(dc_DielMult[7*6+MinusZ][t],d_P3[i]);
			zP=FMUL(dc_DielMult[7*6+PlusZ][t],d_P3[i+pitchXY]);
		
			xP=FADD(xM,xP);
			yP=FADD(yM,yP);
			zP=FADD(zM,zP);
			d_P7[i]=FADD(FADD(xP,yP),FADD(zP,d_P7[i]));
		}
		
		if(iteration%ConvergenceCheck==0)
		{
			double OldTotalEnergy=totalEnergy;
			totalEnergy=0.0;
			EnergyCharge=0.0;
			for(i=0;i<8;i++)
				for(j=0;j<PS.Qnum[i];j++)
			{
				dtmp1=PS.P[i][PS.Qpos[i][j]]*PS.Q[i][j]/PS.Qmult[i][j];
					//ftmp=PS.P[i][PS.Qpos[i][j]]*PS.Q[i][j]/PS.Qmult[i][j];
				EnergyCharge+=dtmp1;
			}
			totalEnergy=EnergyCharge/(fpoh*2.0);
			
			totalChange=totalEnergy-OldTotalEnergy;
			relativeChange=totalChange/totalEnergy;
			
			pnpPrintGroup0("<PoissonIterations Nit=\"%8d\" E=\"%20.16e\" dE=\"%.4e\" rel.E=\"%.4e\" ConvFac=\"%.4e\"/>\n", iteration, totalEnergy, totalChange, relativeChange,ConvFac);
			if(PS.Tolerance!=0.0)
			{
				if(fabs(relativeChange)<=PS.Tolerance)
				{
					printf("Solver has reached the requiered tolerance level\n");
					break;
				}
			}
		}
	}
	m_ContWorld->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////
PoissonSolverOnCuda4Struct* GeneratePoissonSolverOnCuda4Struct(ContWorld* m_ContWorld, int MaxIter,float Rel,float Tol)
{
	PoissonSolverOnCuda4Struct* PS=new PoissonSolverOnCuda4Struct;
	PS->GS.x=m_ContWorld->GridSizeGlobal[0]-2;
	PS->GS.y=m_ContWorld->GridSizeGlobal[1]-2;
	PS->GS.z=m_ContWorld->GridSizeGlobal[2]-2;
	PS->GridScale=m_ContWorld->GridScale;
	PS->Pot=m_ContWorld->Potential;
	PS->MaxIter=MaxIter;
	PS->Rel=Rel;
	PS->om1 = 1.0-PS->Rel;
	PS->om2d6 = PS->Rel/6.0;
	PS->Tol=Tol;
	
	PS->P000=NULL; PS->P100=NULL; PS->P200=NULL; PS->P300=NULL;
	PS->P010=NULL; PS->P110=NULL; PS->P210=NULL; PS->P310=NULL;
	PS->P020=NULL; PS->P120=NULL; PS->P220=NULL; PS->P320=NULL;
	PS->P030=NULL; PS->P130=NULL; PS->P230=NULL; PS->P330=NULL;

	PS->P001=NULL; PS->P101=NULL; PS->P201=NULL; PS->P301=NULL;
	PS->P011=NULL; PS->P111=NULL; PS->P211=NULL; PS->P311=NULL;
	PS->P021=NULL; PS->P121=NULL; PS->P221=NULL; PS->P321=NULL;
	PS->P031=NULL; PS->P131=NULL; PS->P231=NULL; PS->P331=NULL;

	PS->P002=NULL; PS->P102=NULL; PS->P202=NULL; PS->P302=NULL;
	PS->P012=NULL; PS->P112=NULL; PS->P212=NULL; PS->P312=NULL;
	PS->P022=NULL; PS->P122=NULL; PS->P222=NULL; PS->P322=NULL;
	PS->P032=NULL; PS->P132=NULL; PS->P232=NULL; PS->P332=NULL;

	PS->P003=NULL; PS->P103=NULL; PS->P203=NULL; PS->P303=NULL;
	PS->P013=NULL; PS->P113=NULL; PS->P213=NULL; PS->P313=NULL;
	PS->P023=NULL; PS->P123=NULL; PS->P223=NULL; PS->P323=NULL;
	PS->P033=NULL; PS->P133=NULL; PS->P233=NULL; PS->P333=NULL;

	int3 GS=PS->GS;
	PS->ps4GS.x=PS->GS.x/4;
	PS->ps4GS.y=PS->GS.y/4;
	PS->ps4GS.z=PS->GS.z/4;
	int ps4GS_XY=PS->ps4GS.x*PS->ps4GS.y;
	int ps4GS_XYZ=ps4GS_XY*PS->ps4GS.z;
	PS->spltGSWBC.x=PS->ps4GS.x+CUDAXTRAX;
	PS->spltGSWBC.y=PS->ps4GS.y+2;
	PS->spltGSWBC.z=PS->ps4GS.z+2;
	PS->spltGSWBC_X=PS->spltGSWBC.x;
	PS->spltGSWBC_XY=PS->spltGSWBC.x*PS->spltGSWBC.y;
	PS->spltGSWBC_XYZ=PS->spltGSWBC.x*PS->spltGSWBC.y*PS->spltGSWBC.z;
	
	int itmp=PS->spltGSWBC_XYZ;
	PS->P000=new float[itmp]; PS->P100=new float[itmp];
	PS->P200=new float[itmp]; PS->P300=new float[itmp];
	PS->P010=new float[itmp]; PS->P110=new float[itmp];
	PS->P210=new float[itmp]; PS->P310=new float[itmp];
	PS->P020=new float[itmp]; PS->P120=new float[itmp];
	PS->P220=new float[itmp]; PS->P320=new float[itmp];
	PS->P030=new float[itmp]; PS->P130=new float[itmp];
	PS->P230=new float[itmp]; PS->P330=new float[itmp];

	PS->P001=new float[itmp]; PS->P101=new float[itmp];
	PS->P201=new float[itmp]; PS->P301=new float[itmp];
	PS->P011=new float[itmp]; PS->P111=new float[itmp];
	PS->P211=new float[itmp]; PS->P311=new float[itmp];
	PS->P021=new float[itmp]; PS->P121=new float[itmp];
	PS->P221=new float[itmp]; PS->P321=new float[itmp];
	PS->P031=new float[itmp]; PS->P131=new float[itmp];
	PS->P231=new float[itmp]; PS->P331=new float[itmp];

	PS->P002=new float[itmp]; PS->P102=new float[itmp];
	PS->P202=new float[itmp]; PS->P302=new float[itmp];
	PS->P012=new float[itmp]; PS->P112=new float[itmp];
	PS->P212=new float[itmp]; PS->P312=new float[itmp];
	PS->P022=new float[itmp]; PS->P122=new float[itmp];
	PS->P222=new float[itmp]; PS->P322=new float[itmp];
	PS->P032=new float[itmp]; PS->P132=new float[itmp];
	PS->P232=new float[itmp]; PS->P332=new float[itmp];

	PS->P003=new float[itmp]; PS->P103=new float[itmp];
	PS->P203=new float[itmp]; PS->P303=new float[itmp];
	PS->P013=new float[itmp]; PS->P113=new float[itmp];
	PS->P213=new float[itmp]; PS->P313=new float[itmp];
	PS->P023=new float[itmp]; PS->P123=new float[itmp];
	PS->P223=new float[itmp]; PS->P323=new float[itmp];
	PS->P033=new float[itmp]; PS->P133=new float[itmp];
	PS->P233=new float[itmp]; PS->P333=new float[itmp];

	return PS;
}
int FreePoissonSolverOnCuda4Struct(PoissonSolverOnCuda4Struct* PS)
{
	DeleteCArray(PS->P333); DeleteCArray(PS->P233);
	DeleteCArray(PS->P133); DeleteCArray(PS->P033);
	DeleteCArray(PS->P323); DeleteCArray(PS->P223);
	DeleteCArray(PS->P123); DeleteCArray(PS->P023);
	DeleteCArray(PS->P313); DeleteCArray(PS->P213);
	DeleteCArray(PS->P113); DeleteCArray(PS->P013);
	DeleteCArray(PS->P303); DeleteCArray(PS->P203);
	DeleteCArray(PS->P103); DeleteCArray(PS->P003);

	DeleteCArray(PS->P332); DeleteCArray(PS->P232);
	DeleteCArray(PS->P132); DeleteCArray(PS->P032);
	DeleteCArray(PS->P322); DeleteCArray(PS->P222);
	DeleteCArray(PS->P122); DeleteCArray(PS->P022);
	DeleteCArray(PS->P312); DeleteCArray(PS->P212);
	DeleteCArray(PS->P112); DeleteCArray(PS->P012);
	DeleteCArray(PS->P302); DeleteCArray(PS->P202);
	DeleteCArray(PS->P102); DeleteCArray(PS->P002);

	DeleteCArray(PS->P331); DeleteCArray(PS->P231);
	DeleteCArray(PS->P131); DeleteCArray(PS->P031);
	DeleteCArray(PS->P321); DeleteCArray(PS->P221);
	DeleteCArray(PS->P121); DeleteCArray(PS->P021);
	DeleteCArray(PS->P311); DeleteCArray(PS->P211);
	DeleteCArray(PS->P111); DeleteCArray(PS->P011);
	DeleteCArray(PS->P301); DeleteCArray(PS->P201);
	DeleteCArray(PS->P101); DeleteCArray(PS->P001);

	DeleteCArray(PS->P330); DeleteCArray(PS->P230);
	DeleteCArray(PS->P130); DeleteCArray(PS->P030);
	DeleteCArray(PS->P320); DeleteCArray(PS->P220);
	DeleteCArray(PS->P120); DeleteCArray(PS->P020);
	DeleteCArray(PS->P310); DeleteCArray(PS->P210);
	DeleteCArray(PS->P110); DeleteCArray(PS->P010);
	DeleteCArray(PS->P300); DeleteCArray(PS->P200);
	DeleteCArray(PS->P100); DeleteCArray(PS->P000);

	DeleteObjByPnt(PS);
	return 1;
}
#ifndef WITH_CUDA
extern "C" int DoPoissonSolverOnCuda4Float(PoissonSolverOnCuda4Struct* PS)
{
	return 1;
}
extern "C" int DoPoissonSolverOnCuda1Float(PoissonSolverOnCuda1Struct* PS)
{
	return 1;
}
extern "C" int GetCUDADevStat()
{
	return 1;
}
#endif

///////////////////////////////////////////////////////////////////////////////
PoissonSolverOnCuda1Struct* GeneratePoissonSolverOnCuda1Struct(ContWorld* m_ContWorld, int MaxIter,float Rel,float Tol)
{
	PoissonSolverOnCuda1Struct* PS=new PoissonSolverOnCuda1Struct;
	PS->GS.x=m_ContWorld->GridSizeGlobal[0]-2;
	PS->GS.y=m_ContWorld->GridSizeGlobal[1]-2;
	PS->GS.z=m_ContWorld->GridSizeGlobal[2]-2;
	PS->GridScale=m_ContWorld->GridScale;
	PS->Pot=m_ContWorld->Potential;
	PS->MaxIter=MaxIter;
	PS->Rel=Rel;
	PS->om1 = 1.0-PS->Rel;
	PS->om2d6 = PS->Rel/6.0;
	PS->Tol=Tol;
	
	PS->PotCu=NULL;
	
	PS->spltGSWBC.x=PS->GS.x+CUDAXTRAX;
	PS->spltGSWBC.y=PS->GS.y+2;
	PS->spltGSWBC.z=PS->GS.z+2;
	PS->spltGSWBC_X=PS->spltGSWBC.x;
	PS->spltGSWBC_XY=PS->spltGSWBC.x*PS->spltGSWBC.y;
	PS->spltGSWBC_XYZ=PS->spltGSWBC.x*PS->spltGSWBC.y*PS->spltGSWBC.z;
	
	PS->PotCu=new float[PS->spltGSWBC_XYZ];
	int i;
	for(i=0;i<PS->spltGSWBC_XYZ;i++)
		PS->PotCu[i]=0.0f;

	return PS;
}
int FreePoissonSolverOnCuda1Struct(PoissonSolverOnCuda1Struct* PS)
{
	DeleteCArray(PS->PotCu);
	DeleteObjByPnt(PS);
	return 1;
}