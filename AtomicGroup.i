/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester

  This package (LOOS) is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation under version 3 of the License.

  This package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/





%header %{
#include <loos/AtomicGroup.hpp>
#include <loos/sfactories.hpp>
#include <loos/Trajectory.hpp>

#include <loos/Selectors.hpp>
#include <loos/Parser.hpp>
#include <loos/utils.hpp>

#include <loos/Kernel.hpp>
#include <loos/KernelValue.hpp>
#include <loos/KernelActions.hpp>
#include <loos/KernelStack.hpp>
#include <loos/Selectors.hpp>

#include <loos/pdb_remarks.hpp>

#include <sstream>





 
 namespace loos {
   
   
   // This gets translated into Python's StopIteration exception
   struct StopIteration { };
   

   // Iterator class for AtomicGroup
   class AtomicGroupPythonIterator {
   public:
   AtomicGroupPythonIterator(AtomicGroup* p) : _ag(p), _idx(0) { }
     
     pAtom next() throw (loos::StopIteration) {
       if (_idx >= _ag->size())
	 throw(loos::StopIteration());
       return((*_ag)[_idx++]);
     }
   

   private:
     AtomicGroup* _ag;
     uint _idx;
   };
 }


%}


// Translate C++ exception into Python's
%typemap(throws) loos::StopIteration %{
  PyErr_SetNone(PyExc_StopIteration);
  SWIG_fail;
  %}



%wrapper %{
typedef double    greal;
typedef loos::Matrix44<double>   GMatrix;
typedef unsigned int   uint;
typedef unsigned long  ulong;


 %}


namespace loos {


  class AtomicGroupPythonIterator {
  public:
    AtomicGroupPythonIterator(loos::AtomicGroup*);
    loos::pAtom next() throw (loos::StopIteration);
  };
}

%include "loos/AtomicGroup.hpp"

namespace loos {

  %extend AtomicGroup {
    ulong __len__() const {
      return($self->size());
    }


    pAtom __getitem__(const int i) {
      if (i < 0 || static_cast<uint>(i) >= $self->size())
	throw(std::out_of_range("Bad index into AtomicGroup"));
      
      return((*$self)[i]);
    }
    
    void __setitem__(const int i, const pAtom& d) {
      (*$self)[i] = d;
    }

    // Will this leak?
    char* __str__() {
      std::ostringstream oss;
      oss << *$self;
      size_t n = oss.str().size();
      char* buf = new char[n+1];
      strncpy(buf, oss.str().c_str(), n+1);
      return(buf);
      }
    char* __repr__() {
      std::ostringstream oss;
      oss << *$self;
      size_t n = oss.str().size();
      char* buf = new char[n+1];
      strncpy(buf, oss.str().c_str(), n+1);
      return(buf);
    }

    loos::AtomicGroup __copy__() {
      return(loos::AtomicGroup(*$self));
    }

    // We ignore the passed dict...
    loos::AtomicGroup __deepcopy__(PyObject* p) {
      return(loos::AtomicGroup($self->copy()));
    }

    // Bind the AtomicGroup Python iterator to the current group
    loos::AtomicGroupPythonIterator __iter__() {
      return(loos::AtomicGroupPythonIterator($self));
    }
  };

  %rename(__add__)  loos::AtomicGroup::operator+;




};



%template(AtomicGroupVector) std::vector<loos::AtomicGroup>;
