%module loos


%include <std_string.i>
%include <std_vector.i>
%include <boost_shared_ptr.i>

typedef unsigned int    uint;
typedef unsigned long     ulong;

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

