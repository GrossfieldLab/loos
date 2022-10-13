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


%include <std_string.i>

%rename(cpp_splitByMolecule)       loos::AtomicGroup::splitByMolecule;
%rename(cpp_splitByResidue)        loos::AtomicGroup::splitByResidue;
%rename(cpp_splitByUniqueSegid)    loos::AtomicGroup::splitByUniqueSegid;
%rename(cpp_getBondsAGs)           loos::AtomicGroup::getBondsAGs;
%rename(cpp_splitByName)           loos::AtomicGroup::splitByName;



%header %{

#include <map>

#include <AtomicGroup.hpp>
#include <sfactories.hpp>
#include <Trajectory.hpp>

#include <Selectors.hpp>
#include <Parser.hpp>
#include <utils.hpp>

#include <Kernel.hpp>
#include <KernelValue.hpp>
#include <KernelActions.hpp>
#include <KernelStack.hpp>
#include <Selectors.hpp>

#include <pdb_remarks.hpp>

#include <sstream>






 namespace loos {


   // This gets translated into Python's StopIteration exception
   struct StopIteration { };


   // Iterator class for AtomicGroup
   class AtomicGroupPythonIterator {
   public:
   AtomicGroupPythonIterator(AtomicGroup* p) : _ag(p), _idx(0) { }

     pAtom __next__()  {
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
    loos::pAtom __next__();
  };
}


%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* seq, int m, int n)};
%apply (double** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2) {(double** outseq, int* m, int* n)};

%include "AtomicGroup.hpp"

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
    std::string __repr__() const {
      return $self->asString();
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

%pythoncode %{
      def splitByMolecule(self):
          return list(self.cpp_splitByMolecule())

      def splitByResidue(self):
          return list(self.cpp_splitByResidue())

      def splitByUniqueSegid(self):
          return list(self.cpp_splitByUniqueSegid())
      
      def splitByName(self):
          d = {}
          v = self.cpp_splitByName()
          for i in v.keys():
              d[i] = AtomicGroup(v[i])
          return d


      def getBondsAGs(self):
          return list(self.cpp_getBondsAGs())
      
%}
  };

  %rename(__add__)  loos::AtomicGroup::operator+;




};



%template(AtomicGroupVector) std::vector<loos::AtomicGroup>;
%template(AtomicGroupMap)    std::map<std::string, loos::AtomicGroup>;
