/*
  pdb2matlab.cpp

  Takes a PDB and a selection and an optional selection and writes out
  the coordinates to stdout in matlab format...

  The matrix is written in row-major order though, i.e.
  [ X_0 Y_0 Z_0 ;
    X_1 Y_1 Z_1 ;
       ...
    X_i Y_i Z_i ]

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

int main(int argc, char *argv[]) {
  
  if (argc  < 2 || argc > 3) {
    cerr << "Usage: " << argv[0] << " pdb-filename [selection string]" << endl;
    exit(-1);
  }

  PDB pdb(argv[1]);
  AtomicGroup atoms = pdb;

  if (argc > 2)
    atoms = loos::selectAtoms(pdb, argv[2]);

  AtomicGroup::Iterator i(atoms);
  pAtom pa;

  cout << "A = [\n";
  while (pa = i())
    cout << pa->coords().x() << " " << pa->coords().y() << " " << pa->coords().z() << ";\n";
  cout << "];\n";

}
