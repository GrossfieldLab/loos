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




#include "KernelValue.hpp"

namespace loos {

  namespace internal {
  
    int compare(const Value& x, const Value& y) {
      float d;
      int e;

      if (x.type != y.type)
        throw(LOOSError("Comparing values with different types."));

      switch(x.type) {
      case Value::STRING:
        e = (*(x.str) == *(y.str));
        if (!e) {
          if (*(x.str) < *(y.str))
            return(-1);
          else
            return(1);
        }
        return(0);

      case Value::FLOAT:
        d = x.flt - y.flt;
        if (fabs(d) <= FLT_THRESHOLD)
          return(0);
        if (d < 0.0)
          return(-1);
        return(1);
      
      case Value::INT:
        e = x.itg - y.itg;
        return(e);

      case Value::NONE:
      default:
        throw(LOOSError("Invalid comparison"));

      }

    }
  }
}
