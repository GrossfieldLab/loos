/*
  MatrixOrder.hpp

  Matrix ordering policies...
*/

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


#if !defined(LOOS_MATRIX_ORDER_HPP)
#define LOOS_MATRIX_ORDER_HPP

#include <string>
#include <cassert>
#include <stdexcept>

#include <loos/loos_defs.hpp>

namespace loos {

  namespace Math {

    // These are the policy classes for the Matrix class.  They define
    // how the data is actually stored internally, i.e. lower
    // triangular, column major, or row major order...
    //
    // Basically, they store the physical size of the matrix and convert
    // 2D matrix coords to a single index.

  
    //! Class for storing a symmetric triangular matrix
    /**
     * The matrix is lower-triangular.
     */
    class Triangular {
    public:
      explicit Triangular() : m(0), n(0), s(0) { }
      explicit Triangular(const uint y, const uint x) : m(y), n(x), s( (static_cast<ulong>(y)*(y+1))/2 ) {
        assert(y == x && "Cannot have a non-square triangular matrix...  (you know what I mean!)");
      }

      ulong size(void) const { return(s); }

      //! Get the index into the linear array of data
      ulong index(const uint y, const uint x) const {
        ulong b = y;
        ulong a = x;

        if (x > y) {
          b = x;
          a = y;
        }
      
        return( (b * (b + 1)) / 2 + a );
      }

    protected:
      
      //! Reset the [virtual] size of the matrix
      /** Does not currently force a new allocation of data... */
      void setSize(const uint y, const uint x) {
        assert(y == x && "Cannot have a non-square triangular matrix...  (you know what I mean!)");
        m=y; n=x;
        s = (static_cast<ulong>(y)*(y+1))/2;
      }

    
    protected:
      uint m, n;
      ulong s;
    };

    //! Class for storing a matrix in column-major order
    class ColMajor {
    public:
      explicit ColMajor() : m(0), n(0), s(0) { }
      explicit ColMajor(const uint y, const uint x) : m(y), n(x), s(static_cast<ulong>(y)*x) { }

      ulong size(void) const { return(s); }

      //! Get the index into the linear array of data
      ulong index(const uint y, const uint x) const {
        return(static_cast<ulong>(x)*m + y);
      }

    protected:
      //! Reset the [virtual] size of the matrix
      /** Does not currently force a new allocation of data... */
      void setSize(const uint y, const uint x) { m=y; n=x; s = static_cast<ulong>(y)*x; }

    protected:
      uint m, n;
      ulong s;
    };


    //! Class for storing a matrix in row-major order...
    class RowMajor {
    public:
      explicit RowMajor() : m(0), n(0), s(0) { }
      explicit RowMajor(const uint y, const uint x) : m(y), n(x), s(static_cast<ulong>(y)*x) { }

      ulong size(void) const { return(s); }

      //! Get the index into the linear array of data
      ulong index(const uint y, const uint x) const {
        return(static_cast<ulong>(y)*n + x);
      }

    protected:
      //! Reset the [virtual] size of the matrix
      /** Does not currently force a new allocation of data... */
      void setSize(const uint y, const uint x) { m=y; n=x; s=static_cast<ulong>(y)*x; }

    protected:
      uint m, n;
      ulong s;
    };

  }

}

#endif
