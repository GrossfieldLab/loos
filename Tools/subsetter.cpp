/*
  subsetter.cpp

  Subsets a trajectory (stripping out any atoms that don't match the given selection)

  Usage:
    subsetter input-model input-trajectory output-prefix selection-string

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
  string hdr = invocationHeader(argc, argv);

  if (argc != 5) {
    cerr << "Usage- subsetter model trajectory output-prefix selection\n";
    exit(-1);
  }

  AtomicGroup model = createSystem(argv[1]);
  pTraj traj = createTrajectory(argv[2], model);
  string prefix(argv[3]);
  AtomicGroup subset = selectAtoms(model, argv[4]);

  // Now get ready to write the DCD...
  DCDWriter dcdout(prefix + ".dcd");
  dcdout.setTitle(hdr);

  bool first = true;

  while (traj->readFrame()) {
    traj->updateGroupCoords(subset);
    dcdout.writeFrame(subset);

    // Pick off the first frame for the DCD...
    if (first) {
      PDB pdb = PDB::fromAtomicGroup(subset);
      pdb.remarks().add(hdr);
      string out_pdb_name = prefix + ".pdb";
      ofstream ofs(out_pdb_name.c_str());
      ofs << pdb;
      ofs.close();
      first = false;
    }
  }
}
