/*
  Matrix.hpp

  A simple matrix wrapper class.  This is not meant to be a mathematical matrix class!
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



#if !defined(LOOS_MATRIX_HPP)
#define LOOS_MATRIX_HPP

#include <boost/shared_array.hpp>

namespace loos {

  //! Standardize passing of 2-d matrix coords  (row, col)
  struct Duple {
    Duple(const int a, const int b) : j(a), i(b) { }
    Duple() : j(0), i(0) { }

    friend ostream& operator<<(ostream& os, const Duple& d) {
      os << "Duple(" << d.i << "," << d.j << ")";
      return(os);
    }

    int j, i;
  };



  // These are the policy classes for the Matrix class.  They define
  // how the data is actually stored internally, i.e. lower
  // triangular, column major, or row major order...
  // Operator() is where all the magic occurs...  It takes a 2-d
  // matrix coordinate and converts it into a linear index into the
  // array of data...

  //! Class for storing a symmetric triangular matrix
  class Triangular {
  public:
    Triangular(const int y, const int x) : m(y), n(x), s( (y*(y+1))/2 ) {
      assert(y == x && "Cannot have a non-square triangular matrix...  (you know what I mean!)");
    }

    void set(const int y, const int x) {
      assert(y == x && "Cannot have a non-square triangular matrix...  (you know what I mean!)");
      m=y; n=x;
      s = (y*(y+1))/2;
    }

    unsigned long size(void) const { return(s); }

    long operator()(const int y, const int x) const {
      int b = y;
      int a = x;

      if (x > y) {
	b = x;
	a = y;
      }
      
      return( (b * (b + 1)) / 2 + a );
    }

    
  private:
    int m, n;
    long s;
  };

  //! Class for storing a matrix in column-major order
  class ColMajor {
  public:
    ColMajor(const int y, const int x) : m(y), n(x), s(y*x) { }
    void set(const int y, const int x) { m=y; n=x; s = y*x; }
    long size(void) const { return(s); }

    long operator()(const int y, const int x) const { return(x*m + y); }

  private:
    int m, n;
    long s;
  };


  //! Class for storing a matrix in row-major order...
  class RowMajor {
  public:
    RowMajor(const int y, const int x) : m(y), n(x), s(y*x) { }
    void set(const int y, const int x) { m=y; n=x; s=y*x; }
    long size(void) const { return(s); }

    long operator()(const int y, const int x) const { return(y*n + x); }

  private:
    int m, n;
    long s;
  };

  // Reinterpretation functions...
  template<typename T, class C> class Matrix;
  template<typename T> Matrix<T, RowMajor> reinterpretOrder(const Matrix<T, ColMajor>&);
  template<typename T> Matrix<T, ColMajor> reinterpretOrder(const Matrix<T, RowMajor>&);

  // Default order is column major...


  //! Templatized class for intepreting memory as a matrix stored in
  //! different formats...
  /** This class is a wrapper around a block of data allowing the user
   * to access it as a 2D matrix.  The idea is that the data can
   * internally be stored in a row-major, col-major, or triangular
   * format but the user-interface is always the same.  In addition,
   * some capacity to reinterpret the data as col or row major is
   * provided...
   *
   * Internally, the data is managed by a boost shared pointer
   * (boost::shared_array).  Range checking is provided to all
   * accesses.
   *
   * Newly allocate matrices have each element initialized to 0.
   *
   * Note: metadata is currently unused...
   */

  template<typename T, class Policy = ColMajor>
  class Matrix {
  public:

    //! Unitialized matrix
    Matrix() : m(0), n(0), mi(0), pol(0,0), dptr(0), meta("") { }

    //! Wrap an existing block of data with a Matrix.
    Matrix(T* p, const int b, const int a) : m(b), n(a), mi(a*b), pol(b, a), meta("") { dptr = boost::shared_array<T>(p); }

    //! Wrap an already shared block of data...
    Matrix(boost::shared_array<T>& p, const int b, const int a) : m(b), n(a), mi(a*b), pol(b, a), dptr(p), meta("") { }


    //! Create a new block of data for the requested Matrix
    Matrix(const int b, const int a) : m(b), n(a), mi(a*b), pol(b, a), meta("") { allocate(); }

    //! Deep copy of a matrix...
    Matrix<T, Policy> copy(void) const {
      Matrix<T, Policy> result(m, n);
      result.meta = meta;
      long size = pol.size();
      T* p = result.dptr.get();
      T* q = dptr.get();

      for (long i=0; i<size; i++)
	p[i] = q[i];
      
      return(result);
    }

    int rows(void) const { return(m); }
    int cols(void) const { return(n); }

    //! Return a pointer to the internal data...
    T* get(void) { return(dptr.get()); }

    //! Treat the matrix as a 1D array
    T& operator[](const long i) {
      assert(i < mi && "Index out of range in Matrix::operator[]");
      return(dptr[i]);
    }

    //! Return the appropriate element (y-rows, x-cols)
    T& operator()(const int y, const int x) {
      long i = pol(y,x);
      assert(i < mi && "Index out of range in Matrix::operator()");
      return(dptr[i]);
    }

    const T& operator[](const long i) const {
      assert(i < mi && "Index out of range in Matrix::operator[]");
      return(dptr[i]);
    }

    const T& operator()(const int y, const int x) const {
      long i = pol(y,x);
      assert(i < mi && "Index out of range in Matrix::operator()");
      return(dptr[i]);
    }

    void metaData(const string& s) { meta = s; }
    string metaData(void) const { return(meta); }

    //! Deallocate data...
    void free(void) { m = n = 0; dptr.reset(); }

    //! Convert a Col-major to Row-major format
    friend Matrix<T, RowMajor> reinterpretOrder<>(const Matrix<T,ColMajor>&);

    //! Convert a Row-major to Col-major format
    friend Matrix<T, ColMajor> reinterpretOrder<>(const Matrix<T,RowMajor>&);

  private:

    void allocate(void) {
      long size = pol.size();
      dptr.reset();
      dptr = boost::shared_array<T>(new T[size]);
      for (long i=0; i<size; i++)
	dptr[i] = 0;
    }

  private:
    int m, n, mi;
    Policy pol;
    boost::shared_array<T> dptr;
    string meta;
  
  };

  template<typename T>
  Matrix<T, RowMajor> reinterpretOrder(const Matrix<T, ColMajor>& A) {
    Matrix<T, RowMajor> result(A.m, A.n);
    result.dptr = A.dptr;
  
    return(result);
  };

  template<typename T>
  Matrix<T, ColMajor> reinterpretOrder(const Matrix<T, RowMajor>& A) {
    Matrix<T, ColMajor> result(A.m, A.n);
    result.dptr = A.dptr;
  
    return(result);
  };

}

#endif
