/*
  traj2matlab.cpp
*/




/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo
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


#include <loos.hpp>

using namespace std;
using namespace loos;

typedef Math::Matrix<double, Math::RowMajor>   Matrix;



int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  if (argc  != 4) {
    cerr << "Usage: " << argv[0] << " model trajectory selection" << endl;
    exit(-1);
  }

  AtomicGroup model = createSystem(argv[1]);
  pTraj traj = createTrajectory(argv[2], model);
  AtomicGroup subset = selectAtoms(model, argv[3]);

  Matrix A(model.size() * 3, traj->nframes());
  for (uint i = 0; i<traj->nframes(); ++i) {
    traj->readFrame(i);
    traj->updateGroupCoords(subset);
    
    for (int j=0, k=0; j<model.size(); ++j) {
      GCoord c = subset[j]->coords();
      A(k++,i) = c.x();
      A(k++,i) = c.y();
      A(k++,i) = c.z();
    }
  }

  writeAsciiMatrix(cout, A, hdr);
}
