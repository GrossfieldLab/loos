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






#if !defined(MATRIX44_HPP)
#define MATRIX44_HPP


#include <iostream>
#include <stdexcept>
#include <string.h>

#include <Coord.hpp>


using namespace std;


// Forward declaration for matrix-vector multiply
template<class T> Coord<T> operator*(const Matrix44<T>&, const Coord<T>&);


//! Specialized 4x4 Matrix class for handling coordinate transforms.
template<class T>
class Matrix44 {

  T matrix[16];

public:

  //! Create a new identity matrix
  Matrix44() { identity(); }

  //! Create a new matrix with all elements set to v
  Matrix44(const T v) { for (int i = 0; i < 16; i++) matrix[i] = v; }

  //! Zero all elements
  void zero(void) { memset(matrix, 0, 16 * sizeof(T)); }

  //! Identity matrix
  void identity(void) { zero(); matrix[0] = 1; matrix[5] = 1; matrix[10] = 1; matrix[15] = 1; }

  //! Index the matrix element at row j and col i
  T& operator()(const int j, const int i) {
    if (j < 0 || i < 0 || i > 3 || j > 3)
      throw(range_error("Indices into matrix are out of range"));
    return(matrix[j*4+i]);
  }

  //! Index the matrix element at row j and col i
  const T& operator()(const int j, const int i) const {
    if (j < 0 || i < 0 || i > 3 || j > 3)
      throw(range_error("Indices into matrix are out of range"));
    return(matrix[j*4+i]);
  }


  //! Allow access to the linear array of matrix elements
  T& operator[](const int i) {
    if (i < 0 || i > 15)
      throw(range_error("Index into matrix is out of range"));
    return(matrix[i]);
  }

  //! Allow access to the linear array of matrix elements
  const T& operator[](const int i) const {
    if (i < 0 || i > 15)
      throw(range_error("Index into matrix is out of range"));
    return(matrix[i]);
  }

  //! Returns the array pointer
  T* data(void) { return(matrix); }


  //! Addition of two matrices
  Matrix44<T>& operator+=(const Matrix44<T>& rhs) {
    int i;
    for (i=0; i<16; i++)
      matrix[i] += rhs.matrix[i];
    return(*this);
  }

  //! Addition of two matrices
  Matrix44<T> operator+(const Matrix44<T>& rhs) {
    Matrix44<T> res(*this);
    res += rhs;
    return(res);
  }


  //! Addition of a matrix and a constant.
  //! Each element in the matrix is added with the constant...
  //! Relies on the constructor from a constant to handle the case where
  //! you have a matrix + a constant...

  friend Matrix44<T> operator+(const T lhs, const Matrix44<T>& rhs) {
    Matrix44<T> res(rhs);
    res += lhs;
    return(res);
  }

  //! Subtracting matrices
  Matrix44<T>& operator-=(const Matrix44<T>& rhs) {
    int i;
    for (i=0; i<16; i++)
      matrix[i] -= rhs.matrix[i];
    return(*this);
  }

  //! Subtracting matrices
  Matrix44<T> operator-(const Matrix44<T>& rhs) {
    Matrix44<T> res(*this);
    res -= rhs;
    return(res);
  }

  //! Subtraction of a constant from a matrix
  friend Matrix44<T> operator-(const T lhs, const Matrix44<T>& rhs) {
    Matrix44<T> res(rhs);
    res -= lhs;
    return(res);
  }



  //! Friend declaration for matrix-vector multiply...
  friend Coord<T> operator*<>(const Matrix44<T>&, const Coord<T>&);

  //! Matrix-matrix multiply...
  Matrix44<T>& operator*=(const Matrix44<T>& rhs) {
    *this = *this * rhs;
    return(*this);
  }

  //! Matrix-matrix multiply
  Matrix44<T> operator*(const Matrix44<T>& rhs) const {
    Matrix44<T> res;

    res.matrix[0] = matrix[0]*rhs.matrix[0] + matrix[1]*rhs.matrix[4] + matrix[2]*rhs.matrix[8] + matrix[3]*rhs.matrix[12];
    res.matrix[1] = matrix[0]*rhs.matrix[1] + matrix[1]*rhs.matrix[5] + matrix[2]*rhs.matrix[9] + matrix[3]*rhs.matrix[13];
    res.matrix[2] = matrix[0]*rhs.matrix[2] + matrix[1]*rhs.matrix[6] + matrix[2]*rhs.matrix[10] + matrix[3]*rhs.matrix[14];
    res.matrix[3] = matrix[0]*rhs.matrix[3] + matrix[1]*rhs.matrix[7] + matrix[2]*rhs.matrix[11] + matrix[3]*rhs.matrix[15];

    res.matrix[4] = matrix[4]*rhs.matrix[0] + matrix[5]*rhs.matrix[4] + matrix[6]*rhs.matrix[8] + matrix[7]*rhs.matrix[12];
    res.matrix[5] = matrix[4]*rhs.matrix[1] + matrix[5]*rhs.matrix[5] + matrix[6]*rhs.matrix[9] + matrix[7]*rhs.matrix[13];
    res.matrix[6] = matrix[4]*rhs.matrix[2] + matrix[5]*rhs.matrix[6] + matrix[6]*rhs.matrix[10] + matrix[7]*rhs.matrix[14];
    res.matrix[7] = matrix[4]*rhs.matrix[3] + matrix[5]*rhs.matrix[7] + matrix[6]*rhs.matrix[11] + matrix[7]*rhs.matrix[15];

    res.matrix[8] = matrix[8]*rhs.matrix[0] + matrix[9]*rhs.matrix[4] + matrix[10]*rhs.matrix[8] + matrix[11]*rhs.matrix[12];
    res.matrix[9] = matrix[8]*rhs.matrix[1] + matrix[9]*rhs.matrix[5] + matrix[10]*rhs.matrix[9] + matrix[11]*rhs.matrix[13];
    res.matrix[10] = matrix[8]*rhs.matrix[2] + matrix[9]*rhs.matrix[6] + matrix[10]*rhs.matrix[10] + matrix[11]*rhs.matrix[14];
    res.matrix[11] = matrix[8]*rhs.matrix[3] + matrix[9]*rhs.matrix[7] + matrix[10]*rhs.matrix[11] + matrix[11]*rhs.matrix[15];

    res.matrix[12] = matrix[12]*rhs.matrix[0] + matrix[13]*rhs.matrix[4] + matrix[14]*rhs.matrix[8] + matrix[15]*rhs.matrix[12];
    res.matrix[13] = matrix[12]*rhs.matrix[1] + matrix[13]*rhs.matrix[5] + matrix[14]*rhs.matrix[9] + matrix[15]*rhs.matrix[13];
    res.matrix[14] = matrix[12]*rhs.matrix[2] + matrix[13]*rhs.matrix[6] + matrix[14]*rhs.matrix[10] + matrix[15]*rhs.matrix[14];
    res.matrix[15] = matrix[12]*rhs.matrix[3] + matrix[13]*rhs.matrix[7] + matrix[14]*rhs.matrix[11] + matrix[15]*rhs.matrix[15];

    return(res);
  }

  //! Multiplication by a constant.
  /** Each element is multiplied by the same constant
   * This operator [hopefully] prevents auto-casting
   * of the constant to a matrix and then a matrix-matrix
   * multiply, but beware...
   */

  Matrix44<T>& operator*=(const T x) {
    for (int i = 0; i < 16; i++)
      matrix[i] *= x;
    return(*this);
  }

  //! Multiplication by a constant
  Matrix44<T> operator*(const T x) {
    Matrix44<T> res(*this);
    
    res *= x;
    return(res);
  }

  //! Handle the constant * matrix case
  friend Matrix44<T> operator*(const T x, const Matrix44<T>& rhs) {
    Matrix44<T> res(rhs);
    
    res *= x;
    return(res);
  }

  
  //! Output the matrix in pseudo-XML
  friend ostream& operator<<(ostream&os, const Matrix44& m) {
    int i, j, k;

    os << "[";
    for (k=j=0; j<4; j++) {
      for (i=0; i<4; i++)
	os << m[k++] << ((i == 3) ? "" : " ");
      os << ";";
    }
    os << "];";

    return(os);
  }



};


//! Matrix-vector multiply
//! This has to be a friend outside the class for GCC to be happy...
template<class T> Coord<T> operator*(const Matrix44<T>& M, const Coord<T>& v) {
    Coord<T> result;

    result.v[0] = v.v[0]*M.matrix[0] + v.v[1]*M.matrix[1] + v.v[2]*M.matrix[2] + v.v[3]*M.matrix[3];
    result.v[1] = v.v[0]*M.matrix[4] + v.v[1]*M.matrix[5] + v.v[2]*M.matrix[6] + v.v[3]*M.matrix[7];
    result.v[2] = v.v[0]*M.matrix[8] + v.v[1]*M.matrix[9] + v.v[2]*M.matrix[10] + v.v[3]*M.matrix[11];
    result.v[3] = v.v[0]*M.matrix[12] + v.v[1]*M.matrix[13] + v.v[2]*M.matrix[14] + v.v[3]*M.matrix[15];

    return(result);
  }


#endif

