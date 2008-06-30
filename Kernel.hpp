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





#if !defined(KERNEL_HPP)
#define KERNEL_HPP

#include <vector>

#include "Atom.hpp"

#include "KernelValue.hpp"
#include "KernelStack.hpp"
#include "KernelActions.hpp"

using namespace std;

namespace loos {

  //!The Kernel (virtual machine) for compiling and executing user-defined atom selections

  class Kernel {
    vector<Action*> actions;    //! Commands
    ValueStack val_stack;       //! The data stack...

  public:
    
    //! Destroy all of the stored commands...
    ~Kernel() {
      vector<Action*>::iterator i;
      for (i = actions.begin(); i != actions.end(); i++)
	delete (*i);
    }

    //! Add a command, setting the data-stack pointer...
    void push(Action *act) {
      act->setStack(&val_stack);
      actions.push_back(act); 
    }

    void pop(void) { actions.pop_back(); }

    //! Execute the stored commands for a specific atom.
    //! This does not clear the command stack/list
    void execute(pAtom pa = pAtom()) {

      vector<Action*>::iterator i;
      for (i=actions.begin(); i != actions.end(); i++) {
	(*i)->setAtom(pa);
	(*i)->execute();
      }
    }

    
    void clearActions(void) { actions.clear(); }

    ValueStack& stack(void) { return(val_stack); }

    friend ostream& operator<<(ostream& os, const Kernel& k) {
      vector<Action*>::const_iterator i;

      os << "Commands:\n";
      for (i=k.actions.begin(); i != k.actions.end(); i++)
	os << (*i)->name() << endl;
      os << endl;
      return(os);
    }

  };

};


#endif
