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





#if !defined(XFORM_HPP)
#define XFORM_HPP


#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Matrix44.hpp>
#include <Coord.hpp>
#include <loos.hpp>

using namespace std;


typedef Matrix44<greal> GMatrix;

const double PI = 4.0*atan(1.0);



//! Matrix class for handling coordinate transforms...
/** This is based on the OpenGL/RenderMan model of handling geometric
 * transforms.  Coords are expected to be homegenous and the
 * transformation matrix is 4x4.  Rotations are all left-handed.
 *
 * The transform mantains a stack of transformation matrices that you
 * can push and pop as necessary.  You can also load the current
 * transformation with an arbitrary matrix.
 *
 * Transformations are concatenated by post-multiplication.  This means
 * the last declared transformation is the first one applied to an
 * atom's coordinates...  Imagine you've defined a set of
 * transformations in your code: 
 * \verbatim
 *
 *  rotate       ->  M_r
 *  translate    ->  M_t
 *  scale        ->  M_s
 * \endverbatim
 *
 * These are post-multiplied together to create the composite
 * transformation matrix:
 * \verbatim
 *  M = M_r * M_t * M_s
 * \endverbatim
 *
 *  Now, when you transform your coordinate vector, it's just the
 *  matrix-vector multiplication:
 *  \verbatim
 *   v' = Mv = M_r * M_t * M_s * v
 *  \endverbatim
 *
 * So from the perspective of the atom's coordinate frame, we're
 * scaled, then translated, then rotated into the global
 * coordinates... 
 */

class XForm {
  vector<GMatrix> stack;
  bool _unset;

public:
  XForm() : _unset(true) { GMatrix m; stack.push_back(m); }

  //! Push the current matrix onto the stack
  void push(void) { GMatrix M = stack.back(); stack.push_back(M); _unset = false; }
  //! Pop the top matrix off the stack
  void pop(void) {  stack.pop_back(); _unset = false; }
  //! Load a matrix onto the current transform
  void load(const GMatrix& m) { stack.back() = m; _unset = false; }
  //! Concatenate (post-multiply) a matrix with the current transform
  void concat(const GMatrix& m) { stack.back() *= m; _unset = false; }

  //! Premultiply the current transform
  void premult(const GMatrix& m) { GMatrix t = stack.back(); stack.back() = m * t; _unset = false; }

  //! Set the current transform to the identity
  void identity(void) { GMatrix m;  stack.back() = m;  _unset = true; }

  bool unset(void) const { return(_unset); }
  
  //! Translation matrix
  void translate(const greal x, const greal y, const greal z) {
    GMatrix M;

    M(0, 3) = x;
    M(1, 3) = y;
    M(2, 3) = z;
    concat(M);
  }

  //! Translation specified by a GCoord()
  void translate(const GCoord& g) {
    translate(g[0], g[1], g[2]);
  }

  //! Scaling
  void scale(const greal x, const greal y, const greal z) {
    GMatrix M;

    M(0,0) = x;
    M(1,1) = y;
    M(2,2) = z;
    concat(M);
  }

  //! Scaling
  void scale(const GCoord& g) {
    scale(g[0], g[1], g[2]);
  }


  // Angles are in degrees.

  //! Rotate about an arbitrary vector
  //! Angles are specified in degrees.
  void rotate(const GCoord& v, const greal angle) {
    greal theta = PI * angle / 180.0;
    greal c = cos(theta);
    greal s = sin(theta);
    GMatrix M;

    M[0] = v.x() * v.x() * (1.0 - c) + c;
    M[1] = v.x() * v.y() * (1.0 - c) - v.z() * s;
    M[2] = v.x() * v.z() * (1.0 - c) + v.y() * s;

    M[4] = v.x() * v.y() * (1.0 - c) + v.z() * s;
    M[5] = v.y() * v.y() * (1.0 - c) + c;
    M[6] = v.y() * v.z() * (1.0 - c) - v.x() * s;

    M[8] = v.x() * v.z() * (1.0 - c) - v.y() * s;
    M[9] = v.y() * v.z() * (1.0 - c) + v.x() * s;
    M[10] = v.z() * v.z() * (1.0 - c) + c;
    
    concat(M);
  }

  //! Rotate about a specified axis
  //! Axis is given by either 'x', 'y', or 'z'.
  //! Angles are in degrees.
  void rotate(const char axis, const greal angle) {
    switch(axis) {
    case 'x':
    case 'X': rotate(GCoord(1,0,0), angle); break;

    case 'y':
    case 'Y': rotate(GCoord(0,1,0), angle); break;

    case 'z':
    case 'Z': rotate(GCoord(0,0,1), angle); break;

    default:
      throw(logic_error("Invalid axis in XForm::rotate(const char, const greal)"));
    }
  }

  //! Transform a GCoord() with the current transformation
  GCoord transform(const GCoord& v) {
    return(stack.back() * v);
  }

  //! Get the current trasnformation
  // Should we copy or return a ref?
  GMatrix current(void) const {
    GMatrix M = stack.back();
    return(M);
  }


};




#endif
