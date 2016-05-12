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





#if !defined(LOOS_XFORM_HPP)
#define LOOS_XFORM_HPP


#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <cmath>

#include <loos_defs.hpp>
#include <Coord.hpp>
#include <Matrix44.hpp>

namespace loos {

  typedef Matrix44<greal> GMatrix;

  // This was formerly global, but now restricted to this unit...
  namespace {

    const double PI = 4.0*atan(1.0);
    const double very_small = 1e-15;
  }



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
    std::vector<GMatrix> stack;
    bool _unset;

  public:
    XForm() : _unset(true) { GMatrix m; stack.push_back(m); }

    //! Initialize an XForm with an existing matrix.
    explicit XForm(const GMatrix& m) { stack.push_back(m); }

      XForm(const XForm& x) : stack(x.stack), _unset(x._unset) { }

    //! Push the current matrix onto the stack
    void push(void);
    //! Pop the top matrix off the stack
    void pop(void);
    //! Load a matrix onto the current transform
    void load(const GMatrix&);
    //! Concatenate (post-multiply) a matrix with the current transform
    void concat(const GMatrix&);

    //! Premultiply the current transform
    void premult(const GMatrix&);

    //! Set the current transform to the identity
    void identity(void);

    bool unset(void) const;
  
    //! Translation matrix
    void translate(const greal, const greal, const greal);

    //! Translation specified by a GCoord()
    void translate(const GCoord&);

    //! Scaling
    void scale(const greal, const greal, const greal);

    //! Scaling
    void scale(const GCoord&);

    // Angles are in degrees.

    //! Rotate about an arbitrary vector
    //! Angles are specified in degrees.
    void rotate(const GCoord&, const greal);

    //! Rotate about a specified axis
    //! Axis is given by either 'x', 'y', or 'z'.
    //! Angles are in degrees.
    void rotate(const char, const greal);

    //! Transform a GCoord() with the current transformation
    GCoord transform(const GCoord&);

    //! Get the current trasnformation
    // Should we copy or return a ref?
    GMatrix current(void) const;

  };


}

#endif
