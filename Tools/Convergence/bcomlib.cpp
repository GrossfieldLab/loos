/*
  bcomlib.cpp
*/




/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2010, Tod D. Romo
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


#include "bcomlib.hpp"



void Convergence::subtractStructure(loos::RealMatrix& M, const loos::AtomicGroup& model) {
  std::vector<float> avg(model.size() * 3);
  int k = 0;
  for (uint i=0; i<model.size(); ++i) {
    loos::GCoord c = model[i]->coords();
    avg[k++] = c.x();
    avg[k++] = c.y();
    avg[k++] = c.z();
  }

  for (uint i=0; i<M.cols(); ++i)
    for (uint j=0; j<M.rows(); ++j)
      M(j, i) -= avg[j];
}
