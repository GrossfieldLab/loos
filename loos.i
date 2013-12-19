%module loos


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
%include "sfactories.i"
%include "dcdwriter.i"
%include "ensembles.i"
%include "Geometry.i"
%include "TimeSeries.i"
%include "HBondDetector.i"
%include "exceptions.i"
