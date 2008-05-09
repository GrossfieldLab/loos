/*
  XForm.h
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Matrix class for handling coordinate transforms...
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

class XForm {
  vector<GMatrix> stack;

public:
  XForm() { GMatrix m; stack.push_back(m); }

  void push(void) { GMatrix M = stack.back(); stack.push_back(M); }
  void pop(void) {  stack.pop_back(); }
  void load(const GMatrix& m) { stack.back() = m; }
  void concat(const GMatrix& m) { stack.back() *= m; }
  void identity(void) { GMatrix m;  stack.back() = m; }
  
  void translate(const greal x, const greal y, const greal z) {
    GMatrix M;

    M(0, 3) = x;
    M(1, 3) = y;
    M(2, 3) = z;
    concat(M);
  }

  void translate(const GCoord& g) {
    translate(g[0], g[1], g[2]);
  }

  void scale(const greal x, const greal y, const greal z) {
    GMatrix M;

    M(0,0) = x;
    M(1,1) = y;
    M(2,2) = z;
    concat(M);
  }

  void scale(const GCoord& g) {
    scale(g[0], g[1], g[2]);
  }

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

  GCoord transform(const GCoord& v) {
    return(stack.back() * v);
  }

  GMatrix current(void) const {
    GMatrix M = stack.back();
    return(M);
  }


};




#endif
