/*
  paxes.cpp

  Calculates the magnitude of the principal axes (the corresponding
  eigenvalue) over time and writes this out...

  Usage:  paxes mode pdb dcd sel1 sel2 [sel3 ...]
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod Romo
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
#include <boost/format.hpp>

using namespace std;
using namespace loos;


string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\tMagnitudes of the principal axes for a selection over time\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tGiven a selection, the magnitudes of the three principal axes are\n"
    "reported as a function of time.  This gives an idea of the shape of\n"
    "the selection and is a simpler tool to use than molshape.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\tpaxes model.psf trajectory.dcd 'name == \"CA\"' 'resname == \"CAU\"'\n"
    "Reports the time in the first column, followed by the magnitudes of the principal\n"
    "components for all alpha-carbons in the next three columns, followed by the\n"
    "residue named CAU in the following 3 columns.\n"
    "\n"
    "SEE ALSO\n"
    "\tmolshape\n";

  return(msg);
}


int main(int argc, char *argv[]) {
  if (argc < 4) {
    cerr << "Usage - paxes model trajectory sel1 [sel2 ...]\n";
    cerr << fullHelpMessage();
    exit(-1);
  }

  string hdr = invocationHeader(argc, argv);
  AtomicGroup model = createSystem(argv[1]);
  pTraj traj = createTrajectory(argv[2], model);

  cout << "# " << hdr << endl;
  cout << "# frame ";

  vector<AtomicGroup> subsets;
  for (int i=3; i < argc; ++i) {
    AtomicGroup subset = selectAtoms(model, argv[i]);
    subsets.push_back(subset);
    cout << boost::format("a_%d_0 a_%d_1  a_%d_2") % (i-3) % (i-3) % (i-3);
  }
  cout << endl;

  uint t = 0;
  while (traj->readFrame()) {
    traj->updateGroupCoords(model);

    cout << t++ << " ";
    vector<AtomicGroup>::iterator i;
    for (i = subsets.begin(); i != subsets.end(); ++i) {
      vector<GCoord> axes = (*i).principalAxes();
      cout << axes[3][0] << " " << axes[3][1] << " " << axes[3][2] << " ";
    }
    cout << endl;
  }
}


