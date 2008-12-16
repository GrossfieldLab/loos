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





#if !defined(KERNELSTACK_HPP)
#define KERNELSTACK_HPP


#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>


#include <loos_defs.hpp>
#include "KernelValue.hpp"


namespace loos {

  namespace internal {

    class ValueStack {
      std::vector<Value> values;

      void notEmpty(void) const;

    public:
      void push(const Value& val) { values.push_back(val); };

      Value pop(void);

      // Duplicate the top entry...
      void dup(void);

      // Just drop the top entry, i.e. (void)pop()
      void drop(void);

      // Peek at the top value without popping it...
      Value peek(int);

      unsigned int size(void) const { return(values.size()); }

      void clear(void) { values.clear(); }

      friend std::ostream& operator<<(std::ostream&, const ValueStack&);

    };

  }
}


#endif
