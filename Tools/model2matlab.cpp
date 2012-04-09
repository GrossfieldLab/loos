/*
  model2matlab.cpp

  Takes a PDB and a selection and an optional selection and writes out
  the coordinates to stdout in matlab format...

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
  
  cerr << "WARNING- this tool is deprecated and will be removed in a future version of LOOS.\n";
  
  if (argc  != 5) {
    cerr << "Usage: " << argv[0] << " model trajectory selection frame" << endl;
    exit(-1);
  }

  AtomicGroup model = createSystem(argv[1]);
  pTraj traj = createTrajectory(argv[2], model);
  AtomicGroup subset = selectAtoms(model, argv[3]);
  int frame = atoi(argv[4]);

  traj->readFrame(frame);
  traj->updateGroupCoords(subset);
  cout << "A = [\n";
  for (AtomicGroup::iterator i = subset.begin(); i != subset.end(); ++i) {
    GCoord c = (*i)->coords();
    cout << "  " << c.x() << " " << c.y() << " " << c.z() << " ;\n";
  }
  cout << "];\n";

}
