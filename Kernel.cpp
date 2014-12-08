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


#include <Kernel.hpp>
#include <Atom.hpp>

namespace loos {

  Kernel::~Kernel() {
    std::vector<internal::Action*>::iterator i;
    for (i = actions.begin(); i != actions.end(); i++)
      delete (*i);
  }
  
  void Kernel::push(internal::Action *act) {
    act->setStack(&val_stack);
    actions.push_back(act); 
  }
  
  void Kernel::pop(void) {
    if (actions.empty())
      throw(LOOSError("Attempting to pop from an empty action stack"));
    actions.pop_back();
  }
  

  
  void Kernel::execute(pAtom pa) {

    std::vector<internal::Action*>::iterator i;
    for (i=actions.begin(); i != actions.end(); i++) {
      (*i)->setAtom(pa);
      try {
        (*i)->execute();
      }
      catch (LOOSError& e) {
        stack().clear();
        throw(e);
      }
    }
    
  }

    
  void Kernel::clearActions(void) { actions.clear(); }

  internal::ValueStack& Kernel::stack(void) { return(val_stack); }

  std::ostream& operator<<(std::ostream& os, const Kernel& k) {
    std::vector<internal::Action*>::const_iterator i;

    os << "Commands:\n";
    for (i=k.actions.begin(); i != k.actions.end(); i++)
      os << (*i)->name() << std::endl;
    os << std::endl;
    return(os);
  }
  
}

