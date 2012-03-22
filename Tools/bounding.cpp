/*
  bounding.cpp


  Displays the bounding box for a selection from a PDB...
*/


/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo
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

int main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << "Usage: " << argv[0] << " model-filename trajectory selection-string\n"
      "\n"
      "Prints out statistics for the bounding box of the selection over the whole\n"
      "trajectory.  To get the bounding box of a single structure, a PDB may be used as both\n"
      "model and trajectory, i.e. 'bounding foo.pdb foo.dcd all'\n";
    exit(-1);
  }

  AtomicGroup model = createSystem(argv[1]);
  pTraj traj = createTrajectory(argv[2], model);
  AtomicGroup subset = selectAtoms(model, argv[3]);



  GCoord centroid(0,0,0);
  double maxval = numeric_limits<double>::max();
  GCoord max(-maxval, -maxval, -maxval);
  GCoord min(maxval, maxval, maxval);

  while (traj->readFrame()) {
    traj->updateGroupCoords(subset);
    GCoord center(0,0,0);
    for (AtomicGroup::iterator i = subset.begin(); i != subset.end(); ++i) {
      GCoord c = (*i)->coords();
      for (int j = 0; j<3; ++j) {
        if (max[j] < c[j])
          max[j] = c[j];
        if (min[j] > c[j])
          min[j] = c[j];
      }
      center += c;
    }
    center /= subset.size();
    centroid += center;
  }

  centroid /= traj->nframes();

  cout << "Bounds: " << min << " to " << max << endl;
  cout << "Box: " << max - min << endl;
  cout << "Center: " << centroid << endl;
}
