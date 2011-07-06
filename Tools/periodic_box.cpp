/*
  periodic_box.cpp

   Extracts the periodic box information from a trajectory and writes it out to stdout
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2010 Tod D. Romo
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
  
  if (argc != 3) {
    cout << "Usage- " << argv[0] << " model trajectory\n";
    exit(0);
  }

  string hdr = invocationHeader(argc, argv);

  AtomicGroup model = createSystem(argv[1]);
  pTraj traj = createTrajectory(argv[2], model);

  if (! traj->hasPeriodicBox()) {
    cerr << "ERROR- trajectory does not have a periodic box.\n";
    exit(-1);
  }

  cout << "# " << hdr << "\n# t\tX\tY\tZ\n";
  uint t = 0;
  while (traj->readFrame()) {
    GCoord box = traj->periodicBox();
    cout << t++ << "\t" << box.x() << "\t" << box.y() << "\t" << box.z() << endl;
  }
}
