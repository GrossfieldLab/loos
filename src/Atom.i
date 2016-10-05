
%include <std_string.i>
%include <std_vector.i>
%include <boost_shared_ptr.i>

%shared_ptr(loos::Atom)


%{
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>

#include <loos_defs.hpp>
#include "Atom.hpp"

#include <sstream>

%}

namespace loos {
  class Atom;
  typedef boost::shared_ptr<Atom>    pAtom;
}


%include "Atom.hpp"


namespace loos {

  %extend Atom {
    char* __str__() {
      static char buf[1024];
      std::ostringstream oss;
      oss << *$self;
      strncpy(buf, oss.str().c_str(), sizeof(buf));
      return(buf);
    }
    char* __repr__() {
      static char buf[1024];
      std::ostringstream oss;
      oss << *$self;
      strncpy(buf, oss.str().c_str(), sizeof(buf));
      return(buf);
    }

    loos::pAtom __copy__() {
      return(loos::pAtom(new loos::Atom(*$self)));
    }

    // Passed object is ignored...
    loos::pAtom __deepcopy(PyObject* p) {
      return(loos::pAtom(new loos::Atom(*$self)));
    }

    // Only checks equality of principal meta-data, not the atom
    // coordinates...
    bool __eq__(const pAtom& p) {
      return($self->id() == p->id() &&
	     $self->resid() == p->resid() &&
	     $self->resname() == p->resname() &&
	     $self->name() == p->name() &&
	     $self->segid() == p->segid());
    }

    bool __ne__(const pAtom& p) {
      return(! ($self->id() == p->id() &&
                $self->resid() == p->resid() &&
                $self->resname() == p->resname() &&
                $self->name() == p->name() &&
                $self->segid() == p->segid()) );
    }


  };


};

