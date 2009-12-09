/*
  psf-masses

  (c) 2009 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry

  Takes masses from a PSF file and places them into the occupancy field of a PDB.
      
  Notes:

    o Assumes that the atoms are in the same order between the PDB and the PSF
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008,2009 Tod D. Romo
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


using namespace loos;
using namespace std;



int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "Usage- psf-masses model.psf model.pdb >newmodel.pdb\n";
    exit(0);
  }

  string hdr = invocationHeader(argc, argv);

  AtomicGroup source = createSystem(argv[1]);
  AtomicGroup target = createSystem(argv[2]);

  if (source.size() != target.size()) {
    cerr << "ERROR- the files have different number of atoms.\n";
    exit(-1);
  }

  bool flag = false;
  for (int i=0; i<source.size(); ++i) {
    if (source[i]->name() != target[i]->name()) {
      cerr << "ERROR- atom mismatch at position " << i << endl;
      exit(-1);
    }
    if (flag && ! source[i]->checkProperty(Atom::massbit)) {
      flag = false;
      cerr << "WARNING- the PSF does not appear to have masses...using defaults.\n";
    }
    target[i]->occupancy(source[i]->mass());
  }


  PDB pdb = PDB::fromAtomicGroup(target);
  pdb.remarks().add(hdr);
  cout << pdb;
}
