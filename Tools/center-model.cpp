/*
  center-pdb.cpp

  Centers a PDB file at its centroid.
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

  if (argc != 2) {
    cerr << "Usage: center-model structure-file >output.pdb\n";
    exit(-1);
  }

  cerr << "***WARNING: This tool is deprecated and will be removed in a future version of LOOS\n";

  AtomicGroup model = createSystem(argv[1]);
  model.centerAtOrigin();

  PDB pdb = PDB::fromAtomicGroup(model);
  pdb.remarks().add(invocationHeader(argc, argv));
  cout << pdb << endl;
}

