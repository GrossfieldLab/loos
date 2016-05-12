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

#include <XForm.hpp>

namespace loos {

  void XForm::push(void) { GMatrix M = stack.back(); stack.push_back(M); _unset = false; }
  void XForm::pop(void) {  stack.pop_back(); _unset = false; }
  void XForm::load(const GMatrix& m) { stack.back() = m; _unset = false; }
  void XForm::concat(const GMatrix& m) { stack.back() *= m; _unset = false; }

  void XForm::premult(const GMatrix& m) { GMatrix t = stack.back(); stack.back() = m * t; _unset = false; }

  void XForm::identity(void) { GMatrix m;  stack.back() = m;  _unset = true; }

  bool XForm::unset(void) const { return(_unset); }
  
  void XForm::translate(const greal x, const greal y, const greal z) {
    GMatrix M;
    
    M(0, 3) = x;
    M(1, 3) = y;
    M(2, 3) = z;
    concat(M);
  }

  void XForm::translate(const GCoord& g) {
    translate(g[0], g[1], g[2]);
  }

  void XForm::scale(const greal x, const greal y, const greal z) {
    GMatrix M;
    
    M(0,0) = x;
    M(1,1) = y;
    M(2,2) = z;
    concat(M);
  }

  void XForm::scale(const GCoord& g) {
    scale(g[0], g[1], g[2]);
  }


  void XForm::rotate(const GCoord& ov, const greal angle) {
    double l = ov.length();
    if (l < very_small)
      throw(std::invalid_argument("Axis of rotation vector must have non-zero length"));

    GCoord v = ov / ov.length();
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
  
  void XForm::rotate(const char axis, const greal angle) {
    switch(axis) {
    case 'x':
    case 'X': rotate(GCoord(1,0,0), angle); break;
      
    case 'y':
    case 'Y': rotate(GCoord(0,1,0), angle); break;
      
    case 'z':
    case 'Z': rotate(GCoord(0,0,1), angle); break;
      
    default:
      throw(std::logic_error("Invalid axis in XForm::rotate(const char, const greal)"));
    }
  }
  
  GCoord XForm::transform(const GCoord& v) {
    return(stack.back() * v);
  }

  GMatrix XForm::current(void) const {
    GMatrix M = stack.back();
    return(M);
  }

}
