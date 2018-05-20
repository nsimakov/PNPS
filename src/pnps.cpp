/***************************************************************************
 *   Nikolay Simakov                                                       *
 *   nsimakov@andrew.cmu.edu                                               *
 *   Maria Kurnikova Research Group                                        *
 *   http://crete.chem.cmu.edu/                                            *
 *   Carnegie Mellon University, 2005                                      *
 ***************************************************************************/

#ifdef MPI_PARALLEL
  #include <mpi.h>
	#include <unistd.h>
#endif

#if !defined(_MSC_VER) // IGOR FIX to get rid of dependednce on automake config.h
#include "config.h"
#endif

#include "pnpdebug.h"
#include "calcctrl.h"

#include <iostream>
#include <fstream>
#include "tinyxml.h"


// int Run(int argc, char *argv[])
// {
// 
// 	FILE *out=NULL;
// 	if(argc>=4)
// 	{
// 		sprintf(pnptmpdir,"%s\0",argv[3]);
// 	}
// 	else
// 	{
// 		sprintf(pnptmpdir,"/tmp",argv[3]);
// 	}
// 	fprintf(stdout,"pnptmpdir=%s\n",pnptmpdir);		
// 	if(argc>=3)
// 	{
// 		char str[124];
// #ifdef MPI_PARALLEL
// 		sprintf(str,"%s%d\0",argv[2],pnpsapp->GetMyAbsRank());
// #else
// 		sprintf(str,"%s\0",argv[2]);
// #endif
// 		out=fopen(str,"w");
// 		stdout=out;
// 		stderr=out;
// 		setvbuf( stdout, NULL, _IONBF, 0 );
// 	}
// 	if(pnpsapp->GetMyAbsRank()==0)
// 	{
// 		if(argc==1)
// 		{
// 			fprintf(stderr,"ERROR 301: There is not parameters\n");
// 			return EXIT_FAILURE;
// 		}
// 		if(!argv[1])
// 		{
// 			fprintf(stderr,"ERROR 302: Not correct parameters\n");
// 			return EXIT_FAILURE;
// 		}
// 
// 		std::ios::sync_with_stdio();
// 		fprintf(stdout,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\n");
// 		fprintf(stdout,"<Results>\n");
// 		fflush(stdout);
// 	}
// 	
// 	//PNPSapp* app=new PNPSapp();
// 	//app->RunFromFile(argv[1]);
// 	//delete app;app=NULL;
// 	
// 	if(pnpsapp->GetMyAbsRank()==0)
// 	{
// 		fprintf(stdout,"</Results>\n");
// 	}
// 	if(out!=NULL)fclose(out);
// 	return EXIT_SUCCESS;
// }
void PrintHelp()
{
	
}
bool PNPSGetOptBool(const char *OptTag, int argc, char *argv[])
{
	int i;
	for(i=1;i<argc;i++)
	{
		if(strcmp(OptTag,argv[i])==0)
		{
			return true;
		}
	}
	return false;
}
char * PNPSGetOptChar(const char *OptTag, int argc, char *argv[])
{
	int i;
	for(i=1;i<argc-1;i++)
	{
		if(strcmp(OptTag,argv[i])==0)
		{
			if(argv[i+1][0]!='-')
				return argv[i+1];
			else
				return NULL;
		}
	}
	return NULL;
}
int PNPSGetOptPosInt(const char *OptTag, int argc, char *argv[])
{
	int i;
	for(i=1;i<argc-1;i++)
	{
		if(strcmp(OptTag,argv[i])==0)
		{
			if(argv[i+1][0]!='-')
			{
				return atoi(argv[i+1]);
			}
			else
			{
				return -1;
			}
		}
	}
	return -1;
}
TiXmlElement* BcastTiXmlElement(TiXmlElement *Elt)
{
#ifdef MPI_PARALLEL
	int count=0,i;
	char * CStr;
	std::string str;
	
	int MyAbsRank=MPI::COMM_WORLD.Get_rank();
	
	if(MyAbsRank==0)
	{
		str<<*Elt;
		count=str.size();
		
		MPI::COMM_WORLD.Bcast(&count, 1, MPI::INT, 0);
		CStr=new char[count+1];
		for(i=0;i<count;i++)CStr[i]=str[i];
		CStr[count]='\0';
		MPI::COMM_WORLD.Bcast(CStr, count+1, MPI::CHAR, 0);
		delete [] CStr;
		return Elt;
	}
	else
	{
		MPI::COMM_WORLD.Bcast(&count, 1, MPI::INT, 0);
		CStr=new char[count+1];
		MPI::COMM_WORLD.Bcast(CStr, count+1, MPI::CHAR, 0);
		std::istringstream ins(CStr);
		if(Elt==NULL)Elt=new TiXmlElement("DontKnowName");
		ins>>*Elt;
		delete [] CStr;
		return Elt;
	}
#else
	return Elt;
#endif
}
int RunMultiInput(int argc, char *argv[])
{
#ifdef MPI_PARALLEL
	//Read options
	char * MultiInputFile=PNPSGetOptChar("-mi",argc,argv);
	
	PNP_EXIT_FAIL_NULL(MultiInputFile,"no multi input file\n");
	
	
	int MyAbsRank=MPI::COMM_WORLD.Get_rank();
	int TotalNP=MPI::COMM_WORLD.Get_size();
	
	if(MyAbsRank==0)
		InfoPrint("MultiInputFile=%s\n",MultiInputFile);
	
	MPI::COMM_WORLD.Barrier();
	
	//Read parameters from MultiInputFile
	TiXmlElement *RootElt=NULL;
	TiXmlDocument* Doc=NULL;
	if(MyAbsRank==0)
	{
		Doc=new TiXmlDocument(MultiInputFile);
		Doc->LoadFile();
		if(Doc==NULL)
		{
			pnpError("Cannot read file %s\n",MultiInputFile);
			return EXIT_FAILURE;
		}
		RootElt=Doc->RootElement();
		if(RootElt==NULL)
		{
			pnpError("Cannot read file %s\n",MultiInputFile);
			delete Doc;Doc=NULL;
			return EXIT_FAILURE;
		}
		RootElt=BcastTiXmlElement(RootElt);
	}
	else
	{
		RootElt=BcastTiXmlElement(RootElt);
	}
	
	MPI::COMM_WORLD.Barrier();
	//Temp directory
	char * TmpDir=PNPSGetOptChar("-t",argc,argv);
	char pnptmpdir[PNP_MAP_IO_STRING_LENGTH];
	if(TmpDir==NULL)
		sprintf(pnptmpdir,"/tmp");
	else
		sprintf(pnptmpdir,"%s",TmpDir);
	MPI::COMM_WORLD.Barrier();
	
	//Read groups
	int FirstSuperGroupProc=0;
	int GroupNumber=0;
	int countProc=0;
	const TiXmlElement *CldElt;
	CldElt=RootElt->FirstChildElement("Run");
	
	bool GroupLider=false;
	int np=0;
	int ng=0;
	int nppg=0;
	const char * InputFile=NULL;
	const char * WorkDir=NULL;
	const char * OutputFile=NULL;
	while(countProc<TotalNP&&CldElt!=NULL)
	{
		
		InputFile=CldElt->Attribute("Input");
		WorkDir=CldElt->Attribute("WorkDir");
		OutputFile=CldElt->Attribute("Output");
		CldElt->GetIntAttribute("ng",&ng);
		CldElt->GetIntAttribute("nppg",&nppg);
		np=nppg*ng;
		
		if(MyAbsRank>=countProc&&MyAbsRank<countProc+np)
		{
			if(MyAbsRank==countProc)
				GroupLider=true;
			FirstSuperGroupProc=countProc;
			break;
		}
		countProc+=np;
		GroupNumber++;
		CldElt=CldElt->NextSiblingElement("Run");
	}
	
	MPI::COMM_WORLD.Barrier();
	for(int np=0;np<MPI::COMM_WORLD.Get_size();np++)
	{
		if(MPI::COMM_WORLD.Get_rank()==np)
		{
			InfoPrint("MyAbsRank=%d TotalNP=%d\n",MyAbsRank,TotalNP);
			InfoPrint("\tInputFile=%s\n",InputFile);
			InfoPrint("\tWorkDir=%s\n",WorkDir);
			InfoPrint("\tOutputFile=%s\n",OutputFile);
			InfoPrint("\tGroupLider=%d\n",int(GroupLider));
			InfoPrint("\tGroupNumber=%d\n",int(GroupNumber));
			InfoPrint("\tFirstSuperGroupProc=%d\n",int(FirstSuperGroupProc));
		}
		MPI::COMM_WORLD.Barrier();
	};
	//Run groups
	chdir(WorkDir);
	
	
	PNPSApp::InitPNPSApp(np,ng,nppg,pnptmpdir,FirstSuperGroupProc);
	
	PNPSApp::RedirectOutputToFile(OutputFile,false);
	
	const char * cWallTime=PNPSGetOptChar("-walltime",argc,argv);
	if(cWallTime!=NULL)
		pnpsapp->SetWallTimeToRun(cWallTime);
	//print one more time to output file
	InfoPrint("MyAbsRank=%d TotalNP=%d\n",MyAbsRank,TotalNP);
	InfoPrint("\tInputFile=%s\n",InputFile);
	InfoPrint("\tWorkDir=%s\n",WorkDir);
	InfoPrint("\tOutputFile=%s\n",OutputFile);
	InfoPrint("\tGroupLider=%d\n",int(GroupLider));
	InfoPrint("\tGroupNumber=%d\n",int(GroupNumber));
	InfoPrint("\tFirstSuperGroupProc=%d\n",int(FirstSuperGroupProc));
	pnpsapp->RunFromInputFile(InputFile);
	
	
	MPI::COMM_WORLD.Barrier();
	for(int np=0;np<MPI::COMM_WORLD.Get_size();np++)
	{
		if(MPI::COMM_WORLD.Get_rank()==np)
			printf("EVERYTHING DONE %d\n",MPI::COMM_WORLD.Get_rank());
		MPI::COMM_WORLD.Barrier();
	};
	PNPSApp::RedirectOutputToStd();
	PNPSApp::DeletePNPSApp();
	
	if(MyAbsRank==0)
	{
		delete Doc;Doc=NULL;
	}
	else
	{
		delete RootElt;RootElt=NULL;
	}
#endif
	return EXIT_SUCCESS;
}
int Run2(int argc, char *argv[])
{
	//Read options
	
	char * InputFile=PNPSGetOptChar("-i",argc,argv);
	PNP_EXIT_FAIL_NULL(InputFile,"Run2(): no input file\n");
	
	char * OutputFile=PNPSGetOptChar("-o",argc,argv);
	
	char * TmpDir=PNPSGetOptChar("-t",argc,argv);
	
	//Temp directory
	char pnptmpdir[PNP_MAP_IO_STRING_LENGTH];
	if(TmpDir==NULL)
		sprintf(pnptmpdir,"/tmp");
	else
		sprintf(pnptmpdir,"%s",TmpDir);
	
	InfoPrint("InputFile=%s\n",InputFile);
	InfoPrint("pnptmpdir=%s\n",pnptmpdir);
	
	
	//std::ios::sync_with_stdio();
	int np=PNPSGetOptPosInt("-np",argc,argv);
	int ng=PNPSGetOptPosInt("-ng",argc,argv);
	int nppg=PNPSGetOptPosInt("-nppg",argc,argv);
	bool sdmpi=PNPSGetOptBool("-sdmpi",argc,argv);
	bool mdmpi=PNPSGetOptBool("-mdmpi",argc,argv);
	if(sdmpi && mdmpi)
	{
		sdmpi=false;
		mdmpi=false;
	}
	InfoPrint("Befor proccessing: ng=%d; np=%d; nppg=%d;\n",ng,np,nppg);
#ifdef MPI_PARALLEL
	if(np<0)
	{
		np=MPI::COMM_WORLD.Get_size();
	}
	else if(np!=MPI::COMM_WORLD.Get_size())
	{
		pnpError("option -np do not coinside with processors given by MPI\n");
		return EXIT_FAILURE;
	}
	
	if(mdmpi)
		if(ng<0 && nppg<0)
	{
		ng=np;
		nppg=1;
	}
	if(sdmpi)
		if(ng<0 && nppg<0)
	{
		ng=1;
		nppg=np;
	}
	
	if(ng<0 && nppg<0)
	{
		pnpError("option -ng and -nppg needed in MDMPI/SDMPI mixed version\n");
		return EXIT_FAILURE;
	}
	if(ng*nppg!=np)
	{
		pnpError("option -ng, -nppg or -np are not correct(ng*nppg!=np)\n");
		return EXIT_FAILURE;
	}
#else
	if(np<0 && ng<0 && nppg<0)
	{
		ng=1;np=1;nppg=1;
	}
	else if(ng!=1 || np!=1 || nppg!=1)
	{
		ng=1;np=1;nppg=1;
		pnpError("Serial version: -np -ng -nppg are not correct, remove it or use ng=1, np=1 and nppg=1.\n");
		return EXIT_FAILURE;
	}
#endif
	InfoPrint("After proccessing: ng=%d; np=%d; nppg=%d;\n",ng,np,nppg);
	PNPSApp::RedirectOutputToFile(OutputFile,false);
	
	PNPSApp::InitPNPSApp(np,ng,nppg,pnptmpdir);
	const char * cWallTime=PNPSGetOptChar("-walltime",argc,argv);
	if(cWallTime!=NULL)
		pnpsapp->SetWallTimeToRun(cWallTime);
	pnpsapp->RunFromInputFile(InputFile);
	
	
#ifdef MPI_PARALLEL
	MPI::COMM_WORLD.Barrier();
	for(int np=0;np<MPI::COMM_WORLD.Get_size();np++)
	{
		if(MPI::COMM_WORLD.Get_rank()==np)
			printf("EVERYTHING DONE\n");
		MPI::COMM_WORLD.Barrier();
	};
#endif
	PNPSApp::RedirectOutputToStd();
	PNPSApp::DeletePNPSApp();
	return EXIT_SUCCESS;
}

#if !defined(HARLEM_MOD)

int main(int argc, char *argv[])
{
#ifdef MPI_PARALLEL
	MPI::Init(argc, argv);
#endif
	
#ifdef MPI_PARALLEL
	char * MultiInputFile=PNPSGetOptChar("-mi",argc,argv);
	if(MultiInputFile==NULL)
		Run2(argc,argv);
	else
		RunMultiInput(argc,argv);
#else
	Run2(argc,argv);
#endif

#ifdef MPI_PARALLEL
	MPI::Finalize();
#endif
	return EXIT_SUCCESS;
}
#endif