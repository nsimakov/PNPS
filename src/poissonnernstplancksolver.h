//
// C++ Interface: poissonnernstplancksolver
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef POISSONNERNSTPLANCKSOLVER_H
#define POISSONNERNSTPLANCKSOLVER_H

#include <haobject.h>
#include <pnpinterfaces.h>

#ifndef PyObject_HEAD
struct _object;
typedef _object PyObject;
#endif

class ContWorld;
class PoissonSolver;
class NernstPlankSolver;
class VectorField3D;

/**
	@author Nikolay Simakov <nsimakov@andrew.cmu.edu>
*/
class PoissonNernstPlanckSolver : public HaObject, public GenericSolver
{
	public:
		friend class PoissonNernstPlanckMultiGridSolver;
		PoissonNernstPlanckSolver();
		~PoissonNernstPlanckSolver();
		int InitZero();//!<initialize zero values
		int Clear();//!< delete all internal arrays and objects
	public:
		//external variables, parameters for solver
		int MaxIterations;
		/* nernst planck tolerance used to determine pnp convergence,does not work */
		int ConvergenceCheck;
		float tolerance;
		/* dimension alint which to calculate current (x, y, or z) */
		char currentDimension;
		bool verbose;
		int PMFWeightMode;
		std::vector<std::string> PMFWeightModeStr;
		
		//external variables, results
		double Itot,Ipos,Ineg;
		double ItotErr,IposErr,InegErr;
		
		double * positiveCurrentProfile;//! Will be **I in the feture
		double * negativeCurrentProfile;
		
		bool bLimitCurrentCalc;
		float LimitCurrentCalcZ[2];
		bool SaveMemory;
		bool bDouble;
	public:
		//internal variables
		ContWorld *World;
		PoissonSolver *Poisson;
		NernstPlankSolver *NernstPlank;
	public:
		//internal variables
	public:
		//methods
		virtual int SaveXML(TiXmlElement* Elt, HaContext* p_ctxt = NULL );
		virtual int LoadXML(const TiXmlElement* Elt, HaContext* p_ctxt = NULL );
		virtual int LoadParamFromPyDict(PyObject *dict);
		
		int SetContWorld(ContWorld* _world);

		virtual int ShowParameters();//!<show parameters of the solver
    //!show properies can be diffrent from parameters of the solver, valid after InitSolver()
		virtual int ShowProperties(){ return 0;}
    
		virtual int InitSolver();//!<initiate internal arrays
		virtual int Solve();//!< Solve problem
		int SolveSingle();//!< Solve problem
		int SolveDouble();//!< Solve problem
		
		VectorField3D* CalcCartI(int ion);//!< Calculate curent along axis
		VectorField3D* CalcIinout(VectorField3D* CartI);//!< Calculate curent using average from in and out
		VectorField3D* CalcAvrI(VectorField3D* I,int iavr);
	protected:
		int CombineAndPrintCurrents(int iteration);
		int CombineAndPrintCurrentsDouble(int iteration);
};
class PoissonNernstPlanckMultiGridSolver : public HaObject, public GenericSolver
{
	public:
		PoissonNernstPlanckMultiGridSolver();
		~PoissonNernstPlanckMultiGridSolver();
		int InitZero();//!<initialize zero values
		int Clear();//!< delete all internal arrays and objects
	public:
		
		//external variables, results
		double Itot,Ipos,Ineg;
		double ItotErr,IposErr,InegErr;
		
		double * positiveCurrentProfile;//! Will be **I in the feture
		double * negativeCurrentProfile;
		
	public:
		//internal variables
		ContWorld *WorldLev0;
		PoissonNernstPlanckSolver *PNPLev0;
		ContWorld *WorldLev1;
		PoissonNernstPlanckSolver *PNPLev1;
	public:
		//internal variables
	public:
		//methods
		int RunPNPMG(const TiXmlElement* RootElt);
	protected:
		int RunPNPonLev0(const TiXmlElement* Elt);
		int RunPNPonLev1(const TiXmlElement* Elt);
		int ReadPNPnPnNP(const TiXmlElement* Elt);
		
		int DetermineBoarderBetweenLevels();
		int InternalBoarderMinGridLev0[3],InternalBoarderMaxGridLev0[3];
		
	public:
		int RunOnlyPoisson(const TiXmlElement* Elt);
	public:
		virtual int InitSolver();
		virtual int Solve();
		virtual int ShowParameters(){ return 0;}
		virtual int ShowProperties(){ return 0;}
		
};

#ifdef SWIG
%pythoncode %{
def SolvePNPSR(contworld,**kwargs):
	pnps=PoissonNernstPlanckSolver()
	pnps.LoadParamFromPyDict(kwargs);
	pnps.SetContWorld(contworld)
	pnps.InitSolver()
	pnps.Solve()
%}
#endif

#endif
