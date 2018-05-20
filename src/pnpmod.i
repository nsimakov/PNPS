%module pnpmod
%include "exception.i"

%{
#include <cstddef> 
#include <hampi.h>                                                                 
#include "hastl.h"                 
#include "hatypes.h"                                                                    
#include "abstree.h"                     
#include "rapidxml.hpp"                                             
#include "hastring.h"                            
#include "harlemapp.h"                     
#include "hampi.h"           
#include "haxml.h"      
#include "halinalg.h"        
#include "haatom.h"    
#include "hamolecule.h"          
#include "atom_mapping.h" 
#include "hamolset.h"  
#include "moleditor.h"   
#include "haresdb.h"
#include "qc_params.h"
#include "haqchem.h"
#include "hamultipole.h" 
#include "hapseudopot.h"
#include "hagaussian.h"
#include "hazindo.h"
#include "haatbasdb.h"  
#include "halocorb.h"
#include "hacoord.h"
#include "haintcrd.h"
#include "rigidbodycoord.h"
#include "trajanal.h"
#include "haenefunc.h"
#include "haintermol.h"
#include "etcoupl.h"
#include "hamatdb.h"
#include "hatypes.h"
#include "mm_elements.h"
#include "mm_model.h"
#include "mm_params.h"
#include "hamolmech.h" 
#include "mm_traj_anal.h"
#include "mm_driver_amber.h"
#include "mm_driver_tinker.h"
#include "mm_force_field.h"
#include "electrostmod.h"
#include "elmod.h"
#include "nuclacidmod.h"
#include "stmmod.h"
#include "canvas3d.h"
#include "hamolview.h"
#include "haempirical.h"
#include "hasurface.h"
#include "apbsmod.h"
#include "hatests.h"
#include "moleditor.h" 
#include "protonredox.h"
%}

%import ../../../HARLEM/sources/hamolset.i
%import ../../../HALIB/include/halib.i

%include typemaps.i
%include cpointer.i
%include stl.i
%include file.i
%include "carrays.i"
%include "std_vector.i"


%exception { 
    try {
        $action
	} catch (std::runtime_error& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch (...) {
        SWIG_exception(SWIG_RuntimeError, "unknown exception");
    }
}
%{
#include "haxml.h"
#include "pnpmod.h"
#include "mapio.h"
#include "contworld.h"
#include "calcctrl.h"
#include "pmfcalculation.h"
#include "pnpsapp.h"
#include "calcctrl.h"
#include "hasurface.h"
#include "pnpstructs.h"
#include "poissonsolver.h"
#include "poissonboltzmannsolver.h"
#include "poissonnernstplancksolver.h"
#include "pnputil.h"
#include "buildworldni.h"
#include "pbwithljsolver.h"
/*extern "C" {
#include "cps.h"
}*/

#include <stdexcept>

%}


%pointer_class(int,intp);
%pointer_class(double,doublep);
%pointer_functions(double, doublepf);
%array_functions(int, intArray);
%array_functions(bool, boolArray);
%array_functions(float, floatArray);

namespace std {
   %template(vectorf) vector<float>;
};

%include "pnpmod.h"
%include "mapio.h"
%template(HaPyDoubleVectorField3D) HaVectorField3D<double>;

%include "contworld.h"
%include "calcctrl.h"
%include "pmfcalculation.h"
%include "pnpsapp.h"
%include "pnpstructs.h"
%include "poissonsolver.h"
%include "poissonboltzmannsolver.h"
%include "poissonnernstplancksolver.h"
%include "pnputil.h"
%include "buildworldni.h"
%include "pbwithljsolver.h"
/*%pythoncode %{
def addAtoms(builder):
	print "BuildWorldNI::addAtoms()"
%}*/
/*%feature("shadow") BuildWorldNI::addAtoms %{
def addAtoms():
    #do something before
    #$action
    print "BuildWorldNI::addAtoms()"
    #do something after
%}*/

%pythoncode %{
PNPSApp.InitPNPSApp()
%}
/*%include "cps.h"*/

