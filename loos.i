%module loos
%include <exception.i>

%feature("autodoc", "1");

%include <std_string.i>
%include <std_vector.i>
%include <boost_shared_ptr.i>

typedef unsigned int    uint;
typedef unsigned long     ulong;

namespace loos {
  typedef double greal;
  typedef long gint;

  typedef float dcd_real;
  typedef double dcd_double;

  typedef Coord<double> GCoord;
}



# Generic exception wrapper for anything in loos that can look like a list...

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
%include "exceptions.i"
%include "trajwriter.i"
%include "dcdwriter.i"
%include "xtcwriter.i"
%include "sfactories.i"


  %pythoncode %{
from PyTraj import *
	      %}
