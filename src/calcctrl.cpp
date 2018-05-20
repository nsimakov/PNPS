//
// C++ Implementation: pnpsapp
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
#  include "mpi.h"
#endif

#include "calcctrl.h"
#include "contworld.h"
#include "tinyxml.h"
#include "pnpdebug.h"
#include "buildworldni.h"
#include "poissonsolver.h"
#include "sissgausps.h"
#include "poissonboltzmannsolver.h"
#include "pbwithljsolver.h"
#include "poissonnernstplancksolver.h"
#include "nernstplanksolver.h"
#include "pmfcalculation.h"
#include "mapio.h"
#include "pnpconstants.h"
#include "math.h"
#include "pnputil.h"
//#include "dftf.h"

SingleWorldCalcController::SingleWorldCalcController()
{
	InitZero();
}


SingleWorldCalcController::~SingleWorldCalcController()
{
	Clear();
}
int SingleWorldCalcController::InitZero()
{
	m_ContWorld=NULL;
	m_PoissonSolver=NULL;
	m_NernstPlankSolver=NULL;
	m_PoissonNernstPlanckSolver=NULL;
	
	//m_DFTFWorld=NULL;
	
	if(pnpsapp==NULL)PNPSApp::InitPNPSApp();
	return EXIT_SUCCESS;
}
int SingleWorldCalcController::Clear()
{
	DeleteObjByPnt(m_ContWorld);
	DeleteObjByPnt(m_PoissonSolver);
	DeleteObjByPnt(m_PoissonNernstPlanckSolver);
	
	//DeleteObjByPnt(m_DFTFWorld);
	return EXIT_SUCCESS;
}
int SingleWorldCalcController::RunFromCstr(const char *cstr)
{
	TiXmlElement *RootElt=NULL;
	if(pnpsapp->GetMyGroupNumber()==pnpsapp->GetMyGroupLeader())
	{
		TiXmlDocument* Doc=new TiXmlDocument("PNPS");
		Doc->Parse(cstr);
		Doc->Print();
		if(Doc==NULL)
		{
			fprintf(stderr,"Cannot parce document\n");
			return EXIT_FAILURE;
		}
		TiXmlElement *RootElt=Doc->RootElement();
		if(RootElt==NULL)
		{
			fprintf(stderr,"Cannot parce document\n");
			delete Doc;Doc=NULL;
			return EXIT_FAILURE;
		}
		RootElt=pnpsapp->BcastTiXmlElementWithinGroup(RootElt);
		this->Run(RootElt);
		
		delete Doc;Doc=NULL;
	}
	else
	{
		RootElt=pnpsapp->BcastTiXmlElementWithinGroup(RootElt);
		this->Run(RootElt);
		delete RootElt;RootElt=NULL;
	}
	return EXIT_SUCCESS;
}
int SingleWorldCalcController::RunFromFile(const char *FileName)
{
	TiXmlElement *RootElt=NULL;
	if(pnpsapp->GetMyGroupNumber()==pnpsapp->GetMyGroupLeader())
	{
		TiXmlDocument* Doc=new TiXmlDocument(FileName);
		Doc->LoadFile();
		if(Doc==NULL)
		{
			fprintf(stderr,"Cannot read file %s\n",FileName);
			return EXIT_FAILURE;
		}
		TiXmlElement *RootElt=Doc->RootElement();
		if(RootElt==NULL)
		{
			fprintf(stderr,"Cannot read file %s\n",FileName);
			delete Doc;Doc=NULL;
			return EXIT_FAILURE;
		}
		RootElt=pnpsapp->BcastTiXmlElementWithinGroup(RootElt);
		this->Run(RootElt);
		
		delete Doc;Doc=NULL;
	}
	else
	{
		RootElt=pnpsapp->BcastTiXmlElementWithinGroup(RootElt);
		this->Run(RootElt);
		delete RootElt;RootElt=NULL;
	}
	return EXIT_SUCCESS;
}
int SingleWorldCalcController::Run(const TiXmlElement *RootElt)
{
	DefClock0;
	int i,j,k;
	char ctmp0[1024];
	char ctmp1[1024];
	char ctmp2[1024];
	const TiXmlElement *CldElt;
	
	int Status;
	
	if(RootElt==NULL)return EXIT_FAILURE;
	
	CldElt=RootElt->FirstChildElement();
	do
	{
		if(strcmp("PoissonNernstPlanckMultiGridSolver",CldElt->Value())==0||
			strcmp("PNPMGSolver",CldElt->Value())==0)
		{
			PoissonNernstPlanckMultiGridSolver *PNPMGSolver=new PoissonNernstPlanckMultiGridSolver();
			PNPMGSolver->RunPNPMG(CldElt);
			delete PNPMGSolver;
		}
		else if(strcmp("ContWorld",CldElt->Value())==0)
		{
			const TiXmlElement *EltContWorld=CldElt;
			const TiXmlElement *EltContWorldRough=RootElt->FirstChildElement("ContWorldRough");
			const TiXmlElement *EltBuildWorld=RootElt->FirstChildElement("BuildWorldNI");
			const TiXmlElement *EltPoisson=RootElt->FirstChildElement("PoissonSolverRough");
			if(EltContWorld!=NULL && EltContWorldRough!=NULL && EltBuildWorld!=NULL && EltPoisson!=NULL)
				m_ContWorld=CmdContWorldPoissonFocusing(EltContWorld,EltContWorldRough,EltBuildWorld,EltPoisson);
			else
				m_ContWorld=CmdContWorld(CldElt);
		}
		else if(strcmp("BuildWorldNI",CldElt->Value())==0)
		{
			Status=CmdBuildWorld(CldElt,m_ContWorld);
			if(Status==EXIT_FAILURE)return EXIT_FAILURE;
		}
		else if(strcmp("BuildWorldScaled",CldElt->Value())==0)
		{
			Status=CmdBuildWorldScaled(CldElt,m_ContWorld);
			if(Status==EXIT_FAILURE)return EXIT_FAILURE;
		}
		else if(strcmp("MapsIOData",CldElt->Value())==0)
		{
			Status=CmdMapIO(CldElt,m_ContWorld);
			if(Status==EXIT_FAILURE)return EXIT_FAILURE;
		}
		else if(strcmp("PoissonSolver",CldElt->Value())==0)
		{
			m_PoissonSolver=CmdPoissonSolver(CldElt,m_ContWorld);
			PNP_EXIT_FAIL_NULL(m_PoissonSolver,"Fail in Poisson Solver\n");
		}
		else if(strcmp("PoissonSolverOnCuda",CldElt->Value())==0)
		{
			PoissonSolverOnCuda *m_PS_cuda=CmdPoissonSolverOnCuda(CldElt,m_ContWorld);
			delete m_PS_cuda;
		}
		else if(strcmp("PoissonSolverOnCudaDouble",CldElt->Value())==0)
		{
			PoissonSolverOnCudaDouble *m_PS_cuda=CmdPoissonSolverOnCudaDouble(CldElt,m_ContWorld);
			delete m_PS_cuda;
		}
		else if(strcmp("SISSGausPS",CldElt->Value())==0)
		{
			SISSGausPS* m_SISSGausPS=CmdSISSGausPS(CldElt,m_ContWorld);
			PNP_EXIT_FAIL_NULL(m_SISSGausPS,"Fail in SISSGausPS\n");
			delete m_SISSGausPS;
		}
		else if(strcmp("PoissonBoltzmannSolver",CldElt->Value())==0)
		{
			m_PoissonBoltzmannSolver=CmdPoissonBoltzmannSolver(CldElt,m_ContWorld);
			PNP_EXIT_FAIL_NULL(m_PoissonBoltzmannSolver,"Fail in Poisson Boltzmann Solver\n");
		}
		else if(strcmp("PBwithLJSolver",CldElt->Value())==0)
		{
			m_PBwithLJSolver=CmdPBwithLJSolver(CldElt,m_ContWorld);
			PNP_EXIT_FAIL_NULL(m_PBwithLJSolver,"Fail in Poisson Boltzmann with LJ potential Solver\n");
		}
		else if(strcmp("NernstPlankSolver",CldElt->Value())==0)
		{
			m_NernstPlankSolver=CmdNernstPlankSolver(CldElt,m_ContWorld);
			PNP_EXIT_FAIL_NULL(m_NernstPlankSolver,"Fail in Nernst Plank Solver\n");
		}
		else if(strcmp("PoissonNernstPlanckSolver",CldElt->Value())==0)
		{
			m_PoissonNernstPlanckSolver=CmdPoissonNernstPlanckSolver(CldElt,m_ContWorld);
			PNP_EXIT_FAIL_NULL(m_PoissonNernstPlanckSolver,"Fail in Poisson Nernst Plank Solver\n");
			DeleteObjByPnt(m_PoissonNernstPlanckSolver);
		}
		else if(strcmp("PMFProcess",CldElt->Value())==0)
		{
			Status=CmdPMFProcess(CldElt);
			if(Status==EXIT_FAILURE)return EXIT_FAILURE;
		}
		else if(strcmp("AddPotential",CldElt->Value())==0)
		{
			Status=CmdAddPotential(CldElt,m_ContWorld);
			if(Status==EXIT_FAILURE)return EXIT_FAILURE;
		}
		else if(strcmp("AddDiffPotentialToDiffGroup",CldElt->Value())==0)
		{
			Status=CmdAddDiffPotentialToDiffGroup(CldElt,m_ContWorld);
			if(Status==EXIT_FAILURE)return EXIT_FAILURE;
		}
		else if(strcmp("ApplyAsymConc",CldElt->Value())==0)
		{
			Status=CmdApplyAsymConc(CldElt,m_ContWorld);
			if(Status==EXIT_FAILURE)return EXIT_FAILURE;
		}
		else if(strcmp("ApplyAsymConcByScaling",CldElt->Value())==0)
		{
			Status=CmdApplyAsymConcByScaling(CldElt,m_ContWorld);
			if(Status==EXIT_FAILURE)return EXIT_FAILURE;
		}
		else if(strcmp("CheckSystem",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)m_ContWorld->CheckSystem();
		}
		//PNPUtils
		else if(strcmp("SetCatNI",CldElt->Value())==0)
		{
			float NewC;
			if(CldElt->GetFloatAttribute("NewC",&NewC)==EXIT_FAILURE)
			{
				pnpError("No NewC given for SetCatNI\n");
				return EXIT_FAILURE;
			}
			float fpoh= 4*M_PI*m_ContWorld->GridScale;
			float coef=fpoh*COANGS/(m_ContWorld->GridScale*m_ContWorld->GridScale*m_ContWorld->GridScale);
			NewC*=coef;
			Status=CmdSetCatNI(m_ContWorld,NewC);
			if(Status==EXIT_FAILURE)return EXIT_FAILURE;
		}
		else if(strcmp("ReadDynamicChargeZRange",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				std::string filename;
				float Z[2];
				CldElt->GetStdStrAttribute("filename",&filename);
				CldElt->GetArrOfFloatAttribute("Z",Z,2);
				m_ContWorld->ReadDynamicChargeZRange(filename.c_str(),Z[0],Z[1]);
				if(Status==EXIT_FAILURE)return EXIT_FAILURE;
			}
			else
			{
				pnpError("m_ContWorld is not initialized\n");
				return EXIT_FAILURE;
			}
		}
		else if(strcmp("ReadPotentialChargeZRange",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				std::string filename;
				float Z[2];
				CldElt->GetStdStrAttribute("filename",&filename);
				CldElt->GetArrOfFloatAttribute("Z",Z,2);
				m_ContWorld->ReadPotentialChargeZRange(filename.c_str(),Z[0],Z[1]);
				if(Status==EXIT_FAILURE)return EXIT_FAILURE;
			}
			else
			{
				pnpError("m_ContWorld is not initialized\n");
				return EXIT_FAILURE;
			}
		}
		else if(strcmp("PotentialOnBoarderToZero",CldElt->Value())==0)
		{
			Status=CmdPotentialOnBoarderToZero(m_ContWorld);
			if(Status==EXIT_FAILURE)return EXIT_FAILURE;
		}
		else if(strcmp("SetCtoZeroWhereDZero",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				Status=PNPUtil::SetCtoZeroWhereDZero(m_ContWorld);
				if(Status==EXIT_FAILURE)return EXIT_FAILURE;
			}
			else
			{
				pnpError("m_ContWorld is not initialized\n");
				return EXIT_FAILURE;
			}
		}
		else if(strcmp("RemoveQfromNI",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				Status=PNPUtil::RemoveQfromNI(m_ContWorld);
				if(Status==EXIT_FAILURE)return EXIT_FAILURE;
			}
			else
			{
				pnpError("m_ContWorld is not initialized\n");
				return EXIT_FAILURE;
			}
		}
		else if(strcmp("RemoveNodesFromPNPBasedOnNPCriteria",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				float Relaxation;
				if(CldElt->GetFloatAttribute("Relaxation",&Relaxation)==EXIT_FAILURE)Relaxation=1.0;
				Status=PNPUtil::RemoveNodesFromPNPBasedOnNPCriteria(m_ContWorld, Relaxation,5000,-1.0,true);
				if(Status<0)return EXIT_FAILURE;
			}
			else
			{
				pnpError("m_ContWorld is not initialized\n");
				return EXIT_FAILURE;
			}
		}
		else if(strcmp("RemoveDiffusionPointsAtNegativeC",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				Status=m_ContWorld->RemoveDiffusionPointsAtNegativeC();
				if(Status==EXIT_FAILURE)return EXIT_FAILURE;
			}
			else
			{
				pnpError("m_ContWorld is not initialized\n");
				return EXIT_FAILURE;
			}
		}
		else if(strcmp("RemoveDiffusionPoints",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				int n,m_ion;
				int *pnts;
				CldElt->GetIntAttribute("ion",&m_ion);
				CldElt->GetArrOfIntAttributeWithAllocation("pnts",&pnts,&n);
				
				Status=m_ContWorld->RemoveDiffusionPoints(m_ion,pnts,n);
				if(Status==EXIT_FAILURE)return EXIT_FAILURE;
			}
			else
			{
				pnpError("m_ContWorld is not initialized\n");
				return EXIT_FAILURE;
			}
		}
		else if(strcmp("RemoveLargePMFfromPNP",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				float LargePMF;
				CldElt->GetFloatAttribute("LargePMF",&LargePMF);
				Status=PNPUtil::RemoveLargePMFfromPNP(m_ContWorld,LargePMF);
				if(Status==EXIT_FAILURE)return EXIT_FAILURE;
			}
			else
			{
				pnpError("m_ContWorld is not initialized\n");
				return EXIT_FAILURE;
			}
		}
		else if(strcmp("RemoveSmallCfromPNP",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				float SmallC;
				CldElt->GetFloatAttribute("SmallC",&SmallC);
				float fpoh= 4*M_PI*m_ContWorld->GridScale;
				float coef=fpoh*COANGS/(m_ContWorld->GridScale*m_ContWorld->GridScale*m_ContWorld->GridScale);
				SmallC*=coef;
				Status=PNPUtil::RemoveSmallCfromPNP(m_ContWorld,SmallC);
				if(Status==EXIT_FAILURE)return EXIT_FAILURE;
			}
			else
			{
				pnpError("m_ContWorld is not initialized\n");
				return EXIT_FAILURE;
			}
		}
		else if(strcmp("RemoveSmallCfromPNPNearDiffBoarder",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				float SmallC;
				CldElt->GetFloatAttribute("SmallC",&SmallC);
				float fpoh= 4*M_PI*m_ContWorld->GridScale;
				float coef=fpoh*COANGS/(m_ContWorld->GridScale*m_ContWorld->GridScale*m_ContWorld->GridScale);
				SmallC*=coef;
				Status=PNPUtil::RemoveSmallCfromPNPNearDiffBoarder(m_ContWorld,SmallC);
				if(Status==EXIT_FAILURE)return EXIT_FAILURE;
			}
			else
			{
				pnpError("m_ContWorld is not initialized\n");
				return EXIT_FAILURE;
			}
		}
		else if(strcmp("RemoveLargedPMFandLargePMFfromPNP",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				float LargePMF,dPMF;
				CldElt->GetFloatAttribute("LargePMF",&LargePMF);
				CldElt->GetFloatAttribute("dPMF",&dPMF);
				Status=PNPUtil::RemoveLargedPMFandLargePMFfromPNP(m_ContWorld,dPMF,LargePMF);
				if(Status==EXIT_FAILURE)return EXIT_FAILURE;
			}
			else
			{
				pnpError("m_ContWorld is not initialized\n");
				return EXIT_FAILURE;
			}
		}
		else if(strcmp("SetDZeroWithIndex",CldElt->Value())==0)
		{
			Status=CmdSetDZeroWithIndex(m_ContWorld);
			if(Status==EXIT_FAILURE)return EXIT_FAILURE;
		}
		else if(strcmp("SetDNIZeroWithDFlex",CldElt->Value())==0)
		{
			Status=CmdSetDNIZeroWithDFlex(m_ContWorld);
			if(Status==EXIT_FAILURE)return EXIT_FAILURE;
		}
		else if(strcmp("SetInternalCtoZero",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				Status=PNPUtil::SetInternalCtoZero(m_ContWorld);
				if(Status==EXIT_FAILURE)return EXIT_FAILURE;
			}
			else
			{
				pnpError("m_ContWorld is not initialized\n");
				return EXIT_FAILURE;
			}
		}
		else if(strcmp("SetDzeroAtEps",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				int iEps;
				CldElt->GetIntAttribute("iEps",&iEps);
				Status=PNPUtil::SetDzeroAtEps(m_ContWorld,iEps);
				if(Status==EXIT_FAILURE)return EXIT_FAILURE;
			}
			else
			{
				pnpError("m_ContWorld is not initialized\n");
				return EXIT_FAILURE;
			}
		}
		else if(strcmp("ConvertPBLJresultsToDynamicCharge",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				Status=PNPUtil::ConvertPBLJresultsToDynamicCharge(m_ContWorld);
				if(Status==EXIT_FAILURE)return EXIT_FAILURE;
			}
			else
			{
				pnpError("m_ContWorld is not initialized\n");
				return EXIT_FAILURE;
			}
		}
		else if(strcmp("RemoveBadDiffusionPoints",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				Status=m_ContWorld->NIndexing->RemoveBadDiffusionPoints();
				if(Status==EXIT_FAILURE)return EXIT_FAILURE;
			}
			else
			{
				pnpError("m_ContWorld is not initialized\n");
				return EXIT_FAILURE;
			}
		}
		else if(strcmp("RemoveLargedPMFfromPNP",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				float dPMF;
				CldElt->GetFloatAttribute("dPMF",&dPMF);
				Status=PNPUtil::RemoveLargedPMFfromPNP(m_ContWorld,dPMF);
				if(Status==EXIT_FAILURE)return EXIT_FAILURE;
			}
			else
			{
				pnpError("m_ContWorld is not initialized\n");
				return EXIT_FAILURE;
			}
		}
		else if(strcmp("ScaleDiffusionInTheChannel",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				float x0,y0;
				std::vector<float> zRlim;
				std::vector<float> Rlim;
				std::vector<float> zDiffScale;
				std::vector<float> DiffScale;
				int n;
				CldElt->GetFloatAttribute("x0",&x0);
				CldElt->GetFloatAttribute("y0",&y0);
				CldElt->GetTwoVectorFloatElement("Rlim",&zRlim,&Rlim);
				CldElt->GetTwoVectorFloatElement("DiffScale",&zDiffScale,&DiffScale);
				//CldElt->GetTwoArrOfFloatFromText(&Rlim,&DiffScale,&n);
				Status=PNPUtil::ScaleDiffusionInTheChannel(m_ContWorld,  x0, y0,&zRlim,&Rlim,&zDiffScale,&DiffScale);
				if(Status==EXIT_FAILURE)return EXIT_FAILURE;
			}
			else
			{
				pnpError("m_ContWorld is not initialized\n");
				return EXIT_FAILURE;
			}
		}
		else if(strcmp("SaveQstPhi",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				Status=m_ContWorld->SaveQstPhi(CldElt->Attribute("FileName"));
				if(Status==EXIT_FAILURE)return EXIT_FAILURE;
			}
		}
		else if(strcmp("SIPPointsSearch",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				SIPPointsSearch *sipPS=new SIPPointsSearch();
				sipPS->SetContWorld(m_ContWorld);
				sipPS->CmdSIPPointsSearch(CldElt);
				delete sipPS;
			}
		}
		else if(strcmp("ConvertIonStrengthToDynamicCharge",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				m_ContWorld->ConvertIonStrengthToDynamicCharge();
			}
			else
			{
				pnpError("m_ContWorld is not initialize\n");
			}
		}
		else if(strcmp("NPMaskBuilder",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				NPMaskBuilder *m_NPMaskBuilder=new NPMaskBuilder(m_ContWorld);
				m_NPMaskBuilder->CmdNPMaskBuilder(CldElt);
				delete m_NPMaskBuilder;
			}
			else
			{
				pnpError("m_ContWorld is not initialize\n");
			}
		}
		else if(strcmp("IAVCalc",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				IAVCalc *m_IAVCalc=new IAVCalc(m_ContWorld);
				m_IAVCalc->CmdIAVCalc(CldElt);
				delete m_IAVCalc;
			}
			else
			{
				pnpError("m_ContWorld is not initialize\n");
			}
		}
		/*else if(strcmp("DFTFWorld",CldElt->Value())==0)
		{
			m_DFTFWorld=new DFTFWorld();
			m_DFTFWorld->LoadXML(CldElt);
		}
		else if(strcmp("DFTFWorldBuilder",CldElt->Value())==0)
		{
			if(m_DFTFWorld!=NULL)
			{
				DFTFWorldBuilder* blr=new DFTFWorldBuilder();
				blr->LoadXML(CldElt);
				blr->BuildDFTFWorld(m_DFTFWorld);
				delete blr;
			}
			else
			{
				pnpError("m_DFTFWorld is not initialize\n");
			}
		}
		else if(strcmp("DFTFDoSmth",CldElt->Value())==0)
		{
			if(m_DFTFWorld!=NULL)
			{
				DFTFDoSmth* dosmth=new DFTFDoSmth();
				dosmth->LoadXML(CldElt);
				dosmth->DoSmth(m_DFTFWorld);
				delete dosmth;
			}
			else
			{
				pnpError("m_DFTFWorld is not initialize\n");
			}
		}
		else if(strcmp("DFTFSave",CldElt->Value())==0)
		{
			if(m_DFTFWorld!=NULL)
			{
				m_DFTFWorld->SaveMaps(CldElt);
			}
			else
			{
				pnpError("m_DFTFWorld is not initialize\n");
			}
		}*/
		else if(strcmp("CalcRMSDforPotential",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				CmdCalcRMSD(CldElt,m_ContWorld);
			}
			else
			{
				pnpError("m_ContWorld is not initialize\n");
			}
		}
		else if(strcmp("SubtractPotential",CldElt->Value())==0)
		{
			if(m_ContWorld!=NULL)
			{
				CmdSubtractPotential(CldElt,m_ContWorld);
			}
			else
			{
				pnpError("m_ContWorld is not initialize\n");
			}
		}
		else
		{
			
		}
		CldElt=CldElt->NextSiblingElement();
	}
	while(CldElt!=NULL);
	return EXIT_SUCCESS;
}
int SingleWorldCalcController::CmdSetDZeroWithIndex(ContWorld* _ContWorld)
{
	int i,j;
	pnpPrint0("<SingleWorldCalcController::CmdSetDZeroWithIndex>\n");
	for(i=0;i<_ContWorld->NIonsTypes;i++)
	{
		int count=0;
		for(j=0;j<_ContWorld->GS_XYZ;j++)
		{
			if(_ContWorld->NIndexing->GetDiffFloat(i,j)==0.0 && _ContWorld->D[i][j]!=0.0)
			{
				_ContWorld->D[i][j]=0.0;
				count++;
			}
		}
		pnpPrint0("%d point removed for ions %d\n",count,i);
	}
	pnpPrint0("</SingleWorldCalcController::CmdSetDZeroWithIndex>\n");
	return EXIT_SUCCESS;
}
int SingleWorldCalcController::CmdSetDNIZeroWithDFlex(ContWorld* _ContWorld)
{
	int i,j;
	pnpPrint0("<SingleWorldCalcController::CmdSetDNIZeroWithDFlex>\n");
	
	for(i=0;i<_ContWorld->NIonsTypes;i++)
	{
		int count=0;
		for(j=0;j<_ContWorld->GS_XYZ;j++)
		{
			if(_ContWorld->D[i][j]==0.0 && _ContWorld->NIndexing->GetDiffFloat(i,j)!=0.0)
			{
				_ContWorld->NIndexing->SetDiffToZero(i,j);
				count++;
			}
		}
		pnpPrint0("%d point removed for ions %d\n",count,i);
	}
	_ContWorld->NIndexing->CalcDiffBoarder();
	pnpPrint0("</SingleWorldCalcController::CmdSetDNIZeroWithDFlex>\n");
	return EXIT_SUCCESS;
}
int SingleWorldCalcController::CmdSetCatNI(ContWorld* _ContWorld,float NewC)
{
	PNP_EXIT_FAIL_NULL(_ContWorld,"m_ContWorld is not initialized\n");
	PNP_EXIT_FAIL_NULL(_ContWorld->NIndexing,"m_ContWorld->NIndexing is not initialized\n");
	pnpPrint0("<SingleWorldCalcController::CmdSetCatNI>\n");
	int i;
	float fpoh= 4*M_PI*_ContWorld->GridScale;
	float coef= fpoh * COANGS / (_ContWorld->GridScale*_ContWorld->GridScale*_ContWorld->GridScale);
	for(i=0;i<MaxIonTypes;i++)
	{
		if(_ContWorld->NIndexing->C[i]!=0.0)
		{
			pnpPrint0("\tSet C[%d]==%g to %g\n",i,_ContWorld->NIndexing->C[i]/coef,NewC/coef);
			_ContWorld->NIndexing->C[i]=NewC;
		}
	}
	pnpPrint0("</SingleWorldCalcController::CmdSetCatNI>\n");
	return EXIT_SUCCESS;
}
int SingleWorldCalcController::CmdPotentialOnBoarderToZero(ContWorld* _ContWorld)
{
	PNP_EXIT_FAIL_NULL(_ContWorld,"m_ContWorld is not initialized\n");
	PNP_EXIT_FAIL_NULL(_ContWorld->Potential,"m_ContWorld->Potential is not initialized\n");
	pnpPrint0("<SingleWorldCalcController::CmdPotentialOnBoarderToZero>\n");
	int ix,iy,iz,GrdPnt;
	int GS_X=_ContWorld->GridSize[0];
	int GS_Y=_ContWorld->GridSize[1];
	int GS_Z=_ContWorld->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	float *Potential=_ContWorld->Potential;
	//xy
	iz=0;
	for(ix=0;ix<GS_X;ix++)
		for(iy=0;iy<GS_Y;iy++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		Potential[GrdPnt]=0.0;
	}
	iz=GS_Z-2;
	for(ix=0;ix<GS_X;ix++)
		for(iy=0;iy<GS_Y;iy++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		Potential[GrdPnt]=0.0;
	}
	//yz
	ix=0;
	for(iz=0;iz<GS_Z;iz++)
			for(iy=0;iy<GS_Y;iy++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		Potential[GrdPnt]=0.0;
	}
	ix=GS_X-1;
	for(iz=0;iz<GS_Z;iz++)
		for(iy=0;iy<GS_Y;iy++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		Potential[GrdPnt]=0.0;
	}
	iy=0;
	for(iz=0;iz<GS_Z;iz++)
		for(ix=0;ix<GS_X;ix++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		Potential[GrdPnt]=0.0;
	}
	iy=GS_Y-1;
	for(iz=0;iz<GS_Z;iz++)
		for(ix=0;ix<GS_X;ix++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		Potential[GrdPnt]=0.0;
	}
	pnpPrint0("</SingleWorldCalcController::CmdPotentialOnBoarderToZero>\n");
	return EXIT_SUCCESS;
}
ContWorld* SingleWorldCalcController::CmdContWorldPoissonFocusing(const TiXmlElement *EltContWorld,const TiXmlElement *EltContWorldRough,const TiXmlElement *EltBuildWorld,const TiXmlElement *EltPoissonRough)
{
	DbgPrint2("ContWorld* SingleWorldCalcController::CmdContWorldPoissonFocusing\n");
	ContWorld* m_ContWorldRough=CmdContWorld(EltContWorldRough);
	CmdBuildWorld(EltBuildWorld,m_ContWorldRough);
	PoissonSolver *PS=CmdPoissonSolver(EltPoissonRough,m_ContWorldRough);
	delete PS;
	float *tmp[1];
	tmp[0]=m_ContWorldRough->Potential;
	
	VectorField3D * Pot = new VectorField3D(m_ContWorldRough->GridSize, m_ContWorldRough->GridScale, 1, tmp);
	
	m_ContWorldRough->Potential=NULL;
	delete m_ContWorldRough;
	
	ContWorld* m_ContWorldFine=CmdContWorld(EltContWorld);
	
	VectorField3D * PotFine = new VectorField3D(m_ContWorldFine->GridSize, m_ContWorldFine->GridScale, 1);
	//PotFine->FillValue(0.0);
	PotFine->InterpolateFromExt(Pot);
	//Pot->amode=VectorField3D::INTERNAL_ALLOC;
	delete Pot;
	delete [] tmp[0];
	
	m_ContWorldFine->Potential=PotFine->V[0];
	
	PotFine->amode=VectorField3D::EXTERNAL_ALLOC;
	delete PotFine;
	
	return m_ContWorldFine;
}
ContWorld* SingleWorldCalcController::CmdContWorld(const TiXmlElement *Elt)
{
	pnpPrint("<Results%s>\n",Elt->Value());
	
	PNP_EXIT_NULL_ON_NULL(Elt,"CmdContWorld:: TiXmlElement is NULL\n");
	
	ContWorld* _ContWorld=new ContWorld();
	PNP_EXIT_NULL_ON_NULL(_ContWorld,"Can not initiate ContWorld\n");
	
	if(_ContWorld->LoadXML(Elt)==EXIT_FAILURE)
	{
		pnpError("Cannot load ContWorld\n");
		return NULL;
	}
	pnpPrint("</Results%s>\n",Elt->Value());
	return _ContWorld;
}
int SingleWorldCalcController::CmdBuildWorld(const TiXmlElement *Elt,ContWorld* _ContWorld)
{
	DefClock0;
	BuildWorldNI* _BuildWorldNI=NULL;
	fprintf(stdout,"<ResultsBuildWorldNI>\n");
	StartClock0;
	if(_ContWorld!=NULL)
	{
		_BuildWorldNI=new BuildWorldNI();
		_BuildWorldNI->LoadXML(Elt);
		_BuildWorldNI->BuildContWorld(_ContWorld);
		delete _BuildWorldNI;
	}
	else
	{
		fprintf(stderr,"ERROR: There is no grid set yet\n");
	}
	StopClockWMes0("BuildWorldNI");
	fprintf(stdout,"</ResultsBuildWorldNI>\n");
	return EXIT_SUCCESS;
}
int SingleWorldCalcController::CmdBuildWorldScaled(const TiXmlElement *Elt,ContWorld* _ContWorld)
{
	DefClock0;
	BuildWorldScaled* _BuildWorldScaled=NULL;
	fprintf(stdout,"<ResultsBuildWorldScaled>\n");
	StartClock0;
	if(_ContWorld!=NULL)
	{
		_BuildWorldScaled=new BuildWorldScaled();
		_BuildWorldScaled->LoadXML(Elt);
		_BuildWorldScaled->BuildContWorld(_ContWorld);
		delete _BuildWorldScaled;
	}
	else
	{
		fprintf(stderr,"ERROR: There is no grid set yet\n");
	}
	StopClockWMes0("BuildWorldScaled");
	fprintf(stdout,"</ResultsBuildWorldScaled>\n");
	return EXIT_SUCCESS;
}
int SingleWorldCalcController::CmdBuildWorldCmp(const TiXmlElement *Elt,ContWorld* _ContWorld)
{
	DefClock0;
	BuildWorldCmp* _BuildWorldNI=NULL;
	fprintf(stdout,"<ResultsBuildWorldNI>\n");
	StartClock0;
	if(_ContWorld!=NULL)
	{
		_BuildWorldNI=new BuildWorldCmp();
		_BuildWorldNI->LoadXML(Elt);
		_BuildWorldNI->BuildContWorld(_ContWorld);
		delete _BuildWorldNI;
	}
	else
	{
		fprintf(stderr,"ERROR: There is no grid set yet\n");
	}
	StopClockWMes0("BuildWorldNI");
	fprintf(stdout,"</ResultsBuildWorldNI>\n");
	return EXIT_SUCCESS;
}
int SingleWorldCalcController::CmdBuildWorldEu(const TiXmlElement *Elt,ContWorld* _ContWorld)
{
	DefClock0;
	BuildWorldEu* _BuildWorldNI=NULL;
	fprintf(stdout,"<ResultsBuildWorldNI>\n");
	StartClock0;
	if(_ContWorld!=NULL)
	{
		_BuildWorldNI=new BuildWorldEu();
		_BuildWorldNI->LoadXML(Elt);
		_BuildWorldNI->BuildContWorld(_ContWorld);
		delete _BuildWorldNI;
	}
	else
	{
		fprintf(stderr,"ERROR: There is no grid set yet\n");
	}
	StopClockWMes0("BuildWorldNI");
	fprintf(stdout,"</ResultsBuildWorldNI>\n");
	return EXIT_SUCCESS;
}
PoissonSolver* SingleWorldCalcController::CmdPoissonSolver(const TiXmlElement *Elt,ContWorld* _ContWorld)
{
	DefClock0;
	PoissonSolver* _PoissonSolver=NULL;
	fprintf(stdout,"<ResultsPoissonSolver>\n");
	StartClock0;
	if(_ContWorld!=NULL)
	{
		_PoissonSolver=new PoissonSolver();
		_PoissonSolver->LoadXML(Elt);
		_PoissonSolver->SetContWorld(_ContWorld);
		_PoissonSolver->InitSolver();
		_PoissonSolver->Solve();
	}
	else
	{
		fprintf(stderr,"ERROR: There is no grid set yet\n");
	}
	StopClockWMes0("PoissonSolver");
	fprintf(stdout,"</ResultsPoissonSolver>\n");
	return _PoissonSolver;
}
PoissonSolverOnCuda* SingleWorldCalcController::CmdPoissonSolverOnCuda(const TiXmlElement *Elt,ContWorld* _ContWorld)
{
	DefClock0;
	PoissonSolverOnCuda* _PoissonSolver=NULL;
	fprintf(stdout,"<ResultsPoissonSolverOnCuda>\n");
	StartClock0;
	if(_ContWorld!=NULL)
	{
		_PoissonSolver=new PoissonSolverOnCuda();
		_PoissonSolver->LoadXML(Elt);
		_PoissonSolver->SetContWorld(_ContWorld);
		_PoissonSolver->InitSolver();
		_PoissonSolver->Solve();
	}
	else
	{
		fprintf(stderr,"ERROR: There is no grid set yet\n");
	}
	StopClockWMes0("PoissonSolverOnCuda");
	fprintf(stdout,"</ResultsPoissonSolverOnCuda>\n");
	return _PoissonSolver;
}
PoissonSolverOnCudaDouble* SingleWorldCalcController::CmdPoissonSolverOnCudaDouble(const TiXmlElement *Elt,ContWorld* _ContWorld)
{
	DefClock0;
	PoissonSolverOnCudaDouble* _PoissonSolver=NULL;
	fprintf(stdout,"<ResultsPoissonSolverOnCuda>\n");
	StartClock0;
	if(_ContWorld!=NULL)
	{
		_PoissonSolver=new PoissonSolverOnCudaDouble();
		_PoissonSolver->LoadXML(Elt);
		_PoissonSolver->SetContWorld(_ContWorld);
		_PoissonSolver->InitSolver();
		_PoissonSolver->Solve();
	}
	else
	{
		fprintf(stderr,"ERROR: There is no grid set yet\n");
	}
	StopClockWMes0("PoissonSolverOnCuda");
	fprintf(stdout,"</ResultsPoissonSolverOnCuda>\n");
	return _PoissonSolver;
}
int SingleWorldCalcController::CmdCalcRMSD(const TiXmlElement *Elt,ContWorld* _ContWorld)
{
	bool CompWithBC;
	if(Elt->GetBoolAttribute("CompWithBC",&CompWithBC)==EXIT_FAILURE)
		CompWithBC=false;
	return CmdCalcRMSDfromFile(Elt->CStrAttribute("Potential2"),_ContWorld,CompWithBC);
}
int SingleWorldCalcController::CmdCalcRMSDfromFile(const char *filename,ContWorld* _ContWorld,bool CompWithBC)
{
	VectorField3D* VF3D=new VectorField3D(filename);
	
	float *Potential2=VF3D->V[0];
	float *Potential1=_ContWorld->Potential;
	int GrdPnt=1+_ContWorld->GS_X+_ContWorld->GS_XY;
	int ix,iy,iz;
	
	double RMSD=0.0;
	double tmp=fabs(Potential1[GrdPnt]-Potential2[GrdPnt]);
	double tmp2=tmp*tmp;
	double Min=tmp2;
	double Max=tmp2;
	
	if(CompWithBC)
	{
		for(iz=0;iz<_ContWorld->GS_X;iz++)
			for(iy=0;iy<_ContWorld->GS_Y;iy++)
				for(ix=0;ix<_ContWorld->GS_Z;ix++)
		{
			GrdPnt=ix+iy*_ContWorld->GS_X+iz*_ContWorld->GS_XY;
			tmp=Potential1[GrdPnt]-Potential2[GrdPnt];
			tmp2=tmp*tmp;
			RMSD+=tmp2;
			if(Min>tmp2)Min=tmp2;
			if(Max<tmp2)Max=tmp2;
		}
		RMSD/=double(_ContWorld->GS_XYZ);
	}
	else
	{
		for(iz=1;iz<_ContWorld->GS_X-1;iz++)
			for(iy=1;iy<_ContWorld->GS_Y-1;iy++)
				for(ix=1;ix<_ContWorld->GS_Z-1;ix++)
		{
			GrdPnt=ix+iy*_ContWorld->GS_X+iz*_ContWorld->GS_XY;
			tmp=Potential1[GrdPnt]-Potential2[GrdPnt];
			tmp2=tmp*tmp;
			RMSD+=tmp2;
			if(Min>tmp2)Min=tmp2;
			if(Max<tmp2)Max=tmp2;
		}
		RMSD/=double((_ContWorld->GS_X-2)*(_ContWorld->GS_Y-2)*(_ContWorld->GS_Z-2));
	}
	delete VF3D;
	pnpPrint("RMSD with %s is %e. Min is %e. Max is %e\n",filename,sqrt(RMSD),sqrt(Min),sqrt(Max));
	return EXIT_SUCCESS;
}
int SingleWorldCalcController::CmdSubtractPotential(const TiXmlElement *Elt,ContWorld* _ContWorld)
{
	VectorField3D* VF3D=new VectorField3D(Elt->CStrAttribute("Potential2"));
	float *Potential2=VF3D->V[0];
	float *Potential1=_ContWorld->Potential;
	int GrdPnt;
	
	for(GrdPnt=0;GrdPnt<_ContWorld->GS_XYZ;GrdPnt++)
	{
		Potential1[GrdPnt]=Potential1[GrdPnt]-Potential2[GrdPnt];
	}
	return EXIT_SUCCESS;
}
SISSGausPS* SingleWorldCalcController::CmdSISSGausPS(const TiXmlElement *Elt,ContWorld* _ContWorld)
{
	DefClock0;
	SISSGausPS* _PoissonSolver=NULL;
	fprintf(stdout,"<ResultsSISSGausPS>\n");
	StartClock0;
	if(_ContWorld!=NULL)
	{
		_PoissonSolver=new SISSGausPS();
		_PoissonSolver->LoadXML(Elt);
		_PoissonSolver->SetContWorld(_ContWorld);
		_PoissonSolver->InitSolver();
		_PoissonSolver->Solve();
	}
	else
	{
		fprintf(stderr,"ERROR: There is no grid set yet\n");
	}
	StopClockWMes0("PoissonSolver");
	fprintf(stdout,"</ResultsSISSGausPS>\n");
	return _PoissonSolver;
}
PoissonBoltzmannSolver* SingleWorldCalcController::CmdPoissonBoltzmannSolver(const TiXmlElement *Elt,ContWorld* _ContWorld)
{
	DefClock0;
	PoissonBoltzmannSolver* _PoissonBoltzmannSolver=NULL;
	fprintf(stdout,"<ResultsPoissonBoltzmannSolver>\n");
	StartClock0;
	if(_ContWorld!=NULL)
	{
		_PoissonBoltzmannSolver=new PoissonBoltzmannSolver();
		_PoissonBoltzmannSolver->LoadXML(Elt);
		_PoissonBoltzmannSolver->SetContWorld(_ContWorld);
		_PoissonBoltzmannSolver->InitSolver();
		_PoissonBoltzmannSolver->Solve();
	}
	else
	{
		fprintf(stderr,"ERROR: There is no grid set yet\n");
	}
	StopClockWMes0("PoissonBoltzmannSolver");
	fprintf(stdout,"</ResultsPoissonBoltzmannSolver>\n");
	return _PoissonBoltzmannSolver;
}
PBwithLJSolver* SingleWorldCalcController::CmdPBwithLJSolver(const TiXmlElement *Elt,ContWorld* _ContWorld)
{
	DefClock0;
	PBwithLJSolver* _PBwithLJSolver=NULL;
	fprintf(stdout,"<ResultsPBwithLJSolver>\n");
	StartClock0;
	if(_ContWorld!=NULL)
	{
		_PBwithLJSolver=new PBwithLJSolver();
		_PBwithLJSolver->LoadXML(Elt);
		_PBwithLJSolver->SetContWorld(_ContWorld);
		_PBwithLJSolver->InitSolver();
		_PBwithLJSolver->Solve();
	}
	else
	{
		fprintf(stderr,"ERROR: There is no grid set yet\n");
	}
	StopClockWMes0("PBwithLJSolver");
	fprintf(stdout,"</ResultsPBwithLJSolver>\n");
	return _PBwithLJSolver;
}
NernstPlankSolver* SingleWorldCalcController::CmdNernstPlankSolver(const TiXmlElement *Elt,ContWorld* _ContWorld)
{
	DefClock0;
	NernstPlankSolver* _NernstPlankSolver=NULL;
	fprintf(stdout,"<ResultsNernstPlankSolver>\n");
	StartClock0;
	if(_ContWorld!=NULL)
	{
		_NernstPlankSolver=new NernstPlankSolver();
		_NernstPlankSolver->LoadXML(Elt);
		_NernstPlankSolver->SetContWorld(_ContWorld);
		_NernstPlankSolver->InitSolver();
		_NernstPlankSolver->Solve();
	}
	else
	{
		fprintf(stderr,"ERROR: There is no grid set yet\n");
	}
	StopClockWMes0("NernstPlankSolver");
	fprintf(stdout,"</ResultsNernstPlankSolver>\n");
	return _NernstPlankSolver;
}
int SingleWorldCalcController::CmdMapIO(const TiXmlElement *Elt,ContWorld* _ContWorld)
{
	DefClock0;
	fprintf(stdout,"<Results%s>\n",Elt->Value());
	StartClock0;
	PNP_EXIT_FAIL_NULL(_ContWorld,"ContWorld is NULL\n");
	PNP_EXIT_FAIL_NULL(Elt,"TiXmlElement for MapsIOData is NULL\n");
	int status=EXIT_FAILURE;
	if(_ContWorld!=NULL)
	{
		MapsIOData* MapIO=new MapsIOData();
		MapIO->LoadXML(Elt);
		if(MapIO->Mode==MapsIOData::Read)
			status=_ContWorld->ReadMaps(MapIO);
		else if(MapIO->Mode==MapsIOData::Write)
			status=_ContWorld->WriteMaps(MapIO);
		else
		{
			//AutoMode
			//if(Elt->
		}
		delete MapIO;
	}
	else
	{
		fprintf(stderr,"ERROR: There is no grid set yet\n");
	}
	PNP_EXIT_ON_FAIL_MES(status,"Can't read or write\n");
	StopClockWMes0("MapsIOData");
	fprintf(stdout,"</Results%s>\n",Elt->Value());
	return EXIT_SUCCESS;
}
int SingleWorldCalcController::CmdAddPotential(const TiXmlElement *Elt,ContWorld* _ContWorld)
{
	if(Elt==NULL)
		{pnpPrint0("ERROR: TiXmlElement is not present.");return EXIT_FAILURE;}
	if(_ContWorld==NULL)
		{pnpPrint0("ERROR: ContWorld not init.");return EXIT_FAILURE;}
	float z[2],phi[2];
	Elt->GetArrOfFloatAttribute("Z",z,2);
	Elt->GetArrOfFloatAttribute("Phi",phi,2);
	return _ContWorld->AddPotential(z[0],z[1],phi[0],phi[1]);
}
int SingleWorldCalcController::CmdAddDiffPotentialToDiffGroup(const TiXmlElement *Elt,ContWorld* _ContWorld)
{
	if(Elt==NULL)
	{pnpPrint0("ERROR: TiXmlElement is not present.");return EXIT_FAILURE;}
	if(_ContWorld==NULL)
	{pnpPrint0("ERROR: ContWorld not init.");return EXIT_FAILURE;}
	float z[2];
	Elt->GetArrOfFloatAttribute("Z",z,2);
	int ndphi=Elt->GetNumOfElmOfArrAttribute("Phi");
	if(ndphi%2==1)
	{
		pnpError("odd number of potential, need to be Ngroups*2\n");
		ndphi=(ndphi-1)/2;
	}
	else
	{
		ndphi=ndphi/2;
	}
	if(ndphi!=pnpsapp->GetNumberOfGroups())
		pnpError("Number of potential diffrences is not equil to number of groups\n");
	
	DbgPrint2("ndphi=%d\n",ndphi);
	
	float *phi=new float[2*ndphi];
	Elt->GetArrOfFloatAttribute("Phi",phi,2*ndphi);
	int result=_ContWorld->AddPotential(z[0],z[1],phi[2*pnpsapp->GetMyGroupNumber()],phi[2*pnpsapp->GetMyGroupNumber()+1]);
	delete [] phi;	
	return result;
}
int SingleWorldCalcController::CmdApplyAsymConcByScaling(const TiXmlElement *Elt,ContWorld* _ContWorld)
{
	if(Elt==NULL)
	{pnpPrint0("ERROR: TiXmlElement is not present.");return EXIT_FAILURE;}
	if(_ContWorld==NULL)
	{pnpPrint0("ERROR: ContWorld not init.");return EXIT_FAILURE;}
	
	float z;
	float *Cz0;
	float *Cz1;
	int i;
	Elt->GetFloatAttribute("Z",&z);
	Elt->GetArrOfFloatAttributeWithAllocation("Cz0",&Cz0,&i);
	if(i!=_ContWorld->NIonsTypes)
		pnpError("ApplyAsymConc number of inputs for Cz0 is not coinside with number of ions type in ContWorld\n");
	Elt->GetArrOfFloatAttributeWithAllocation("Cz1",&Cz1,&i);
	if(i!=_ContWorld->NIonsTypes)
		pnpError("ApplyAsymConc number of inputs for Cz0 is not coinside with number of ions type in ContWorld\n");
	int status=_ContWorld->ApplyAsymConcByScaling(z,Cz0,Cz1);
	delete [] Cz0;
	delete [] Cz1;
	return status;
}
int SingleWorldCalcController::CmdApplyAsymConc(const TiXmlElement *Elt,ContWorld* _ContWorld)
{
	if(Elt==NULL)
	{pnpPrint0("ERROR: TiXmlElement is not present.");return EXIT_FAILURE;}
	if(_ContWorld==NULL)
	{pnpPrint0("ERROR: ContWorld not init.");return EXIT_FAILURE;}
	
	float z;
	float *Cz0;
	float *Cz1;
	int i;
	Elt->GetFloatAttribute("Z",&z);
	Elt->GetArrOfFloatAttributeWithAllocation("Cz0",&Cz0,&i);
	if(i!=_ContWorld->NIonsTypes)
		pnpError("ApplyAsymConc number of inputs for Cz0 is not coinside with number of ions type in ContWorld\n");
	Elt->GetArrOfFloatAttributeWithAllocation("Cz1",&Cz1,&i);
	if(i!=_ContWorld->NIonsTypes)
		pnpError("ApplyAsymConc number of inputs for Cz0 is not coinside with number of ions type in ContWorld\n");
	int status=_ContWorld->ApplyAsymConcByScaling(z,Cz0,Cz1);
	delete [] Cz0;
	delete [] Cz1;
	return status;
}
int SingleWorldCalcController::CmdPMFProcess(const TiXmlElement *Elt)
{
	if(m_ContWorld!=NULL)
	{
		DefClock0;
		fprintf(stdout,"<ResultsPMFProcess>\n");
		StartClock0;
		PMFProcessing* pmfproc=new PMFProcessing();
		pmfproc->SetContWorld(m_ContWorld);
		pmfproc->CmdPMFProcess(Elt);
		delete pmfproc;
		StopClockWMes0("PMFProcess");
		fprintf(stdout,"</ResultsPMFPMFProcess>\n");
		return EXIT_SUCCESS;
	}
	else
	{
		fprintf(stderr,"m_ContWorld is not initialize\n");
		return EXIT_FAILURE;
	}
}
PoissonNernstPlanckSolver*  SingleWorldCalcController::CmdPoissonNernstPlanckSolver(const TiXmlElement *Elt,ContWorld* _ContWorld)
{
	PoissonNernstPlanckSolver*  _PoissonNernstPlanckSolver=NULL;
	DefClock0;
	PoissonSolver* _PoissonSolver=NULL;
	fprintf(stdout,"<ResultsPoissonNernstPlanckSolver>\n");
	StartClock0;
	if(_ContWorld!=NULL)
	{
		_PoissonNernstPlanckSolver=new PoissonNernstPlanckSolver();
		_PoissonNernstPlanckSolver->LoadXML(Elt);
		_PoissonNernstPlanckSolver->SetContWorld(_ContWorld);
		_PoissonNernstPlanckSolver->InitSolver();
		_PoissonNernstPlanckSolver->Solve();
		
		std::string filename;
		if(Elt->GetStdStrAttribute("SaveJ0",&filename)==EXIT_SUCCESS)
		{
			pnpPrint("Save j0 to %s\n",filename.c_str());
			float qrev=1.0/_ContWorld->IonsQ[0];
			pnpPrint("qrev=%f\n",qrev);
			bool avr;
			if(Elt->GetBoolAttribute("AvrJ",&avr)!=EXIT_SUCCESS)
				avr=false;
			VectorField3D* VF=_PoissonNernstPlanckSolver->CalcCartI(0);
			//VF->MultiplyBy(qrev);
			if(avr)
			{
				VectorField3D* VF2=_PoissonNernstPlanckSolver->CalcAvrI(VF,0);
				delete VF;
				VF=VF2;
			}
			VF->WriteToFile(filename.c_str());
			delete VF;
		}
		if(Elt->GetStdStrAttribute("SaveJ1",&filename)==EXIT_SUCCESS)
		{
			pnpPrint("Save j1 to %s\n",filename.c_str());
			float qrev=1.0/_ContWorld->IonsQ[1];
			pnpPrint("qrev=%f\n",qrev);
			bool avr;
			if(Elt->GetBoolAttribute("AvrJ",&avr)!=EXIT_SUCCESS)
				avr=false;
			VectorField3D* VF=_PoissonNernstPlanckSolver->CalcCartI(1);
			//VF->MultiplyBy(qrev);
			if(avr)
			{
				VectorField3D* VF2=_PoissonNernstPlanckSolver->CalcAvrI(VF,0);
				delete VF;
				VF=VF2;
			}
			VF->WriteToFile(filename.c_str());
			delete VF;
		}
	}
	else
	{
		fprintf(stderr,"ERROR: There is no grid set yet\n");
	}
	StopClockWMes0("PoissonNernstPlanckSolver");
	fprintf(stdout,"</ResultsPoissonNernstPlanckSolver>\n");
	return _PoissonNernstPlanckSolver;
}
