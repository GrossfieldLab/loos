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



string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\tPrint statistics for bounding box of a selection\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tPrint out statistics for the bounding box of a selection over the whole\n"
    "trajectory.  To get the bounding box of a single structure, a PDB may be used\n"
    "as both model and trajectory.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\tbounding model.psf simulation.dcd 'name == \"CA\"'\n"
    "Bounding box for all alpha-carbons\n"
    "\n"
    "\tbounding model.pdb model.pdb 'name == \"CA\"'\n"
    "Bounding box for a single structure.\n"
    "\n"
    "NOTES\n"
    "\tThe bounding box of a model (no trajectory) ONLY works for PDB files\n"
    "\n";

  return(msg);
}


int main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << "Usage: " << argv[0] << " model-filename trajectory selection-string\n";
    cerr << fullHelpMessage();
    exit(-1);
  }

  AtomicGroup model = createSystem(argv[1]);
  pTraj traj = createTrajectory(argv[2], model);
  AtomicGroup subset = selectAtoms(model, argv[3]);



  GCoord centroid(0,0,0);
  double maxval = numeric_limits<double>::max();
  GCoord max(-maxval, -maxval, -maxval);
  GCoord min(maxval, maxval, maxval);
  vector<GCoord> boxes;
  GCoord avgbox;

  while (traj->readFrame()) {
    traj->updateGroupCoords(subset);

    GCoord center(0,0,0);
    GCoord submin(maxval, maxval, maxval);
    GCoord submax(-maxval, -maxval, -maxval);

    for (AtomicGroup::iterator i = subset.begin(); i != subset.end(); ++i) {
      GCoord c = (*i)->coords();
      for (int j = 0; j<3; ++j) {
        if (submax[j] < c[j])
          submax[j] = c[j];
        if (submin[j] > c[j])
          submin[j] = c[j];
      }
      center += c;
    }
    center /= subset.size();
    centroid += center;

    GCoord box = submax - submin;
    boxes.push_back(box);
    avgbox += box;

    for (uint i=0; i<3; ++i) {
      if (submax[i] > max[i])
	max[i] = submax[i];
      if (submin[i] < min[i])
	min[i] = submin[i];
    }
  }

  centroid /= traj->nframes();
  avgbox /= traj->nframes();

  GCoord boxdev;
  for (vector<GCoord>::const_iterator i = boxes.begin(); i != boxes.end(); ++i) {
    GCoord d = *i - avgbox;
    for (uint j=0; j<3; ++j)
      d[j] *= d[j];
    boxdev += d;
  }
  for (uint j=0; j<3; ++j)
    boxdev[j] = sqrt(boxdev[j]/(boxes.size()-1));

  cout << "Bounds: " << min << " to " << max << endl;
  cout << "Average Box: " << avgbox << endl;
  cout << "Stddev Box: " << boxdev << endl;
  cout << "Center: " << centroid << endl;
}
