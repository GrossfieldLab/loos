%module loos

%include <exception.i>

%feature("autodoc", "1");

%{
#define SWIG_FILE_WITH_INIT
%}

%include <std_string.i>
%include <std_vector.i>
%include <boost_shared_ptr.i>


%include "numpy.i"
%init %{
  import_array();
%}


typedef unsigned int    uint;
typedef unsigned long     ulong;

namespace loos {
  typedef double greal;
  typedef long gint;

  typedef float dcd_real;
  typedef double dcd_double;

  typedef Coord<double> GCoord;
}


/*
 Generic exception wrapper for anything in loos that can look like a list...
*/

%exception __getitem__
{
  try {
    $action
      }
  catch (std::out_of_range& e) {
    PyErr_SetString(PyExc_IndexError, const_cast<char*>(e.what()));
    return(NULL);
  }
}

%template(DoubleVector)         std::vector<double>;
%template(DoubleVectorMatrix)   std::vector< std::vector<double> >;
%template(IntVector)            std::vector<int>;
%template(UIntVector)           std::vector<uint>;
%template(StringVector)         std::vector<std::string>;


%include "exceptions.i"

%include "catch_it.i"

%include "Coord.i"
%include "Atom.i"
%include "Matrix44.i"
%include "pdb_remarks.i"
%include "XForm.i"
%include "AtomicGroup.i"
%include "Trajectory.i"
%include "utils.i"
%include "cryst.i"
%include "pdb.i"
%include "ensembles.i"
%include "Geometry.i"
%include "TimeSeries.i"
%include "HBondDetector.i"
%include "trajwriter.i"
%include "dcdwriter.i"
%include "xtcwriter.i"
%include "sfactories.i"
%include "alignment.i"
%include "gro.i"
%include "utils_structural.i"
%include "utils_random.i"
%include "Weights.i"
%include "RnaSuite.i"
%include "FormFactorSet.i"
