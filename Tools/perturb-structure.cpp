/*
  perturb-structure.cpp

  Apply a random perturbation to a structure (random directions, fixed magnitude)
*/



/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2010, Tod D. Romo
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
    cout << "Usage- perturb_structure magnitude [selection|all] structure-file >output.pdb\n";
    exit(-1);
  }

  
  string hdr = invocationHeader(argc, argv);

  int k = 1;
  double magnitude = strtod(argv[k++], 0);
  string selection(argv[k++]);
  AtomicGroup model = createSystem(argv[k++]);

  AtomicGroup subset = selectAtoms(model, selection);

  subset.perturbCoords(magnitude);
  PDB pdb = PDB::fromAtomicGroup(model);
  pdb.remarks().add(hdr);

  cout << pdb;

}
