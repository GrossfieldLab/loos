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

#include <boost/shared_ptr.hpp>

namespace loos {

  // These are the policy classes for the Matrix class.  They define
  // how the data is actually stored internally, i.e. lower
  // triangular, column major, or row major order...

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

    struct iterator {
      iterator(const int y, const int x) : n((y*(y+1))/2), i(0) { }
      long next(void) {
	if (i >= n)
	  return(-1);
	return(i++);
      }

      long n, i;
    };

    
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


    struct iterator {
      iterator(const int y, const int x) : m(y), n(x), a(0), b(0), eod(false) { }
      long next(void) {
	if (eod)
	  return(-1);

	long i = a*m + b;
	if (++a >= n) {
	  a = 0;
	  if (++b >= m)
	    eod = true;
	}
	return(i);
      }

      int m, n, a, b;
      bool eod;
    };
	  

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


    struct iterator {
      iterator(const int y, const int x) : m(y), n(x), a(0), b(0), eod(false) { }
      long next(void) {
	if (eod)
	  return(-1);

	long i = b*n + a;
	if (++a >= n) {
	  a = 0;
	  if (++b >= m)
	    eod = true;
	}
	return(i);
      }

      int m, n, a, b;
      bool eod;
    };


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
   * (boost::shared_ptr).  Range checking is provided to all
   * accesses.
   *
   * Newly allocate matrices have each element initialized to 0.
   */

  template<typename T, class Policy = ColMajor>
  class Matrix {
  public:

    //! Wrap an existing block of data with a Matrix.
    Matrix(T* p, const int b, const int a) : m(b), n(a), mi(a*b), index(b, a) { dptr = boost::shared_ptr<T>(p); }
    //! Create a new block of data for the requested Matrix
    Matrix(const int b, const int a) : m(b), n(a), mi(a*b), index(b, a) { allocate(); }
    
    //! Deep copy of a matrix...
    Matrix<T> copy(void) const {
      Matrix<T> result(m, n);
      long size = index.size();
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
      T* p = dptr.get();
      assert(i < mi && "Index out of range in Matrix::operator[]");
      return(p[i]);
    }

    //! Return the appropriate element (y-rows, x-cols)
    T& operator()(const int y, const int x) {
      T* p = dptr.get();
      long i = index(y,x);
      assert(i < mi && "Index out of range in Matrix::operator()");
      return(p[i]);
    }

    const T& operator[](const long i) const {
      T* p = dptr.get();
      assert(i < mi && "Index out of range in Matrix::operator[]");
      return(p[i]);
    }

    const T& operator()(const int y, const int x) const {
      T* p = dptr.get();
      long i = index(y,x);
      assert(i < mi && "Index out of range in Matrix::operator()");
      return(p[i]);
    }

    //! Deallocate data...
    void free(void) { m = n = 0; dptr.reset(); }

    //! Convert a Col-major to Row-major format
    friend Matrix<T, RowMajor> reinterpretOrder<>(const Matrix<T,ColMajor>&);

    //! Convert a Row-major to Col-major format
    friend Matrix<T, ColMajor> reinterpretOrder<>(const Matrix<T,RowMajor>&);

    struct iterator {
      iterator(const Matrix<T, Policy>& M) : miter(M.m, M.n) { }
      long next(void) { return(miter.next()); }

      typename Policy::iterator miter;
    };


  private:

    void allocate(void) {
      long size = index.size();
      dptr.reset();
      dptr = boost::shared_ptr<T>(new T[size]);
      T* p = dptr.get();
      for (long i=0; i<size; i++)
	p[i] = 0;
    }

  private:
    int m, n, mi;
    Policy index;
    boost::shared_ptr<T> dptr;
  
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

};

#endif
