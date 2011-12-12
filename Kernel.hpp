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





#if !defined(LOOS_KERNEL_HPP)
#define LOOS_KERNEL_HPP

#include <vector>

#include "KernelValue.hpp"
#include "KernelStack.hpp"
#include "KernelActions.hpp"

namespace loos {

  //!The Kernel (virtual machine) for compiling and executing user-defined atom selections

  class Kernel {
    std::vector<internal::Action*> actions;    //! Commands
    internal::ValueStack val_stack;       //! The data stack...

  public:
    
    //! Destroy all of the stored commands...
    ~Kernel();

    //! Add a command, setting the data-stack pointer...
    void push(internal::Action*);

    void pop(void);

    //! Execute the stored commands for a specific atom.
    //! This does not clear the command stack/list
    void execute(pAtom pa = pAtom());
    
    void clearActions(void);

    internal::ValueStack& stack(void);

    friend std::ostream& operator<<(std::ostream&, const Kernel&);
  };
};


#endif
