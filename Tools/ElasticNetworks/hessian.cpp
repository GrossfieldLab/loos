/*
  hessian

  (c) 2009,2010 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry

*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009 Tod D. Romo
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


#include "hessian.hpp"


using namespace loos;


DoubleMatrix DistanceCutoff::blockImpl(const uint i, const uint j) {

  DoubleMatrix B(3,3);
  GCoord u = nodes[i]->coords();
  GCoord v = nodes[j]->coords();
  GCoord d = v - u;

  double s = d.length2();
  if (s <= radius) {

    for (int j=0; j<3; ++j)
      for (int i=0; i<3; ++i)
        B(i,j) = d[i]*d[j] / s;
  }

  return(B);
}


DoubleMatrix DistanceWeight::blockImpl(const uint i, const uint j) {

  DoubleMatrix B(3,3);
  GCoord u = nodes[i]->coords();
  GCoord v = nodes[j]->coords();
  GCoord d = v - u;

  double s = d.length();
  s = pow(s, power);
  for (int j=0; j<3; ++j)
    for (int i=0; i<3; ++i)
      B(i,j) = d[i]*d[j] * s;

  return(B);
}



DoubleMatrix hessian(SuperBlock* blockMethod) {
  
  uint n = blockMethod->size();
  DoubleMatrix H(3*n,3*n);

  for (uint i=1; i<n; ++i) {
    for (uint j=0; j<i; ++j) {
      DoubleMatrix B = blockMethod->block(i, j);
      for (uint x = 0; x<3; ++x)
        for (uint y = 0; y<3; ++y) {
          H(i*3 + y, j*3 + x) = -B(y, x);
          H(j*3 + x, i*3 + y) = -B(x ,y);
        }
    }
  }

  // Now handle the diagonal...
  for (uint i=0; i<n; ++i) {
    DoubleMatrix B(3,3);
    for (uint j=0; j<n; ++j) {
      if (j == i)
        continue;
      
      for (uint x=0; x<3; ++x)
        for (uint y=0; y<3; ++y)
          B(y,x) += H(j*3 + y, i*3 + x);
    }

    for (uint x=0; x<3; ++x)
      for (uint y=0; y<3; ++y)
        H(i*3 + y, i*3 + x) = -B(y,x);
  }

  return(H);
}

