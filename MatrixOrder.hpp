/*
  MatrixOrder.hpp

  Matrix ordering policies...
*/

/*
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester
*/



#if !defined(LAB_MATRIX_ORDER_HPP)
#define LAB_MATRIX_ORDER_HPP

#include <string>
#include <cassert>
#include <stdexcept>

namespace lab {


  // These are the policy classes for the Matrix class.  They define
  // how the data is actually stored internally, i.e. lower
  // triangular, column major, or row major order...
  //
  // Basically, they store the physical size of the matrix and convert
  // 2D matrix coords to a single index.

  
  //! Class for storing a symmetric triangular matrix
  class Triangular {
  public:
    explicit Triangular() : m(0), n(0), s(0) { }
    explicit Triangular(const uint y, const uint x) : m(y), n(x), s( (y*(y+1))/2 ) {
      assert(y == x && "Cannot have a non-square triangular matrix...  (you know what I mean!)");
    }

    ulong size(void) const { return(s); }

    ulong index(const uint y, const uint x) const {
      if (x >= n || y >= m)
        throw(std::out_of_range("Matrix indices out of range"));

      uint b = y;
      uint a = x;

      if (x > y) {
	b = x;
	a = y;
      }
      
      return( (b * (b + 1)) / 2 + a );
    }

  protected:
    void setSize(const uint y, const uint x) {
      assert(y == x && "Cannot have a non-square triangular matrix...  (you know what I mean!)");
      m=y; n=x;
      s = (y*(y+1))/2;
    }

    
  protected:
    uint m, n;
    ulong s;
  };

  //! Class for storing a matrix in column-major order
  class ColMajor {
  public:
    explicit ColMajor() : m(0), n(0), s(0) { }
    explicit ColMajor(const uint y, const uint x) : m(y), n(x), s(y*x) { }

    ulong size(void) const { return(s); }

    ulong index(const uint y, const uint x) const {
      if (x >= n || y >= m)
        throw(std::out_of_range("Matrix indices out of range"));
      return(x*m + y);
    }

  protected:
    void setSize(const uint y, const uint x) { m=y; n=x; s = y*x; }

  protected:
    uint m, n;
    ulong s;
  };


  //! Class for storing a matrix in row-major order...
  class RowMajor {
  public:
    explicit RowMajor() : m(0), n(0), s(0) { }
    explicit RowMajor(const uint y, const uint x) : m(y), n(x), s(y*x) { }

    ulong size(void) const { return(s); }

    ulong index(const uint y, const uint x) const {
      if (x >= n || y >= m)
        throw(std::out_of_range("Matrix indices out of range"));
      return(y*n + x);
    }

  protected:
    void setSize(const uint y, const uint x) { m=y; n=x; s=y*x; }

  protected:
    uint m, n;
    ulong s;
  };

}



#endif
