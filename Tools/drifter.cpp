/*
  drifter.cpp

  Calculates the average centroid for a selection, then writes out the
  distance between the centroid of each frame's selection and the
  average
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
#include <string>
#include <vector>

using namespace std;
using namespace loos;


GCoord averageCentroid(AtomicGroup& model, pTraj traj) {
  GCoord avg(0,0,0);
  uint n = traj->nframes();
  for (uint i=0; i<n; ++i) {
    traj->readFrame(i);
    traj->updateGroupCoords(model);
    GCoord c = model.centroid();
    avg += c;
  }

  avg /= n;
  return(avg);
}


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  if (argc != 4) {
    cerr << "Usage: drifter selection model trajectory\n";
    exit(-1);
  }

  AtomicGroup model = createSystem(argv[2]);
  pTraj traj = createTrajectory(argv[3], model);

  AtomicGroup subset = selectAtoms(model, argv[1]);

  GCoord avg_center = averageCentroid(subset, traj);

  cout << "# " << hdr << endl;
  cout << "# t d\n";
  uint n = traj->nframes();
  for (uint i=0; i<n; ++i) {
    traj->readFrame(i);
    traj->updateGroupCoords(subset);
    GCoord c = subset.centroid();
    cout << i << " " << avg_center.distance(c) << endl;
  }
}
