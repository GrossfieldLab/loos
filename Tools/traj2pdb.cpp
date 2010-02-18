/*
  traj2dcd


  Converts a LOOS-supported format to a series of PDB files

  Usage:

    traj2pdb model-file trajectory-file pdb-corename

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

#include <iostream>
#include <sstream>
#include <loos.hpp>
#include <boost/format.hpp>

using namespace std;
using namespace loos;


int main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << "Usage - traj2pdb model trajectory pdb-corename\n";
    exit(-1);
  }

  AtomicGroup model = createSystem(argv[1]);
  pTraj traj = createTrajectory(argv[2], model);
  string pdb_core = string(argv[3]);

  uint n = traj->nframes();
  PDB pdb = PDB::fromAtomicGroup(model);
  pdb.remarks().add(invocationHeader(argc, argv));

  cout << boost::format("There are %d atoms and %d frames.\n") % model.size() % n;

  cout << "Processing - ";
  cout.flush();

  for (uint i=0; i<n; ++i) {
    if (i % 250 == 0) {
      cout << '.';
      cout.flush();
    }
    traj->readFrame(i);
    traj->updateGroupCoords(model);

    stringstream s;
    s << pdb_core << "_" << i << ".pdb";
    ofstream pdbout(s.str().c_str());
    if (pdbout.fail())
        {
        cerr << "Error writing file " 
             << s.str()
             << ".  Exiting"
             << endl;
        exit(-1);
        }

    pdbout << pdb;
    pdbout.close();
  }
  
  cout << " done\n";
}

