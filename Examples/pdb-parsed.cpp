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




#include <pdb.hpp>
#include <Parser.hpp>
#include <Selectors.hpp>


int main(int argc, char *argv[]) {
  
  PDB p(argv[1]);

  cout << "Read in " << p.size() << " atoms from " << argv[1] << endl;

  Parser parsed(argv[2]);
  KernelSelector parsed_sel(parsed.kernel());

  cout << "*** Virtual Machine Command STACK ***\n" << parsed.kernel() << endl;
  AtomicGroup  usel = p.select(parsed_sel);
  
  cout << "There are " << usel.size() << " atoms in the selection.\n";
  cout << "The max radius is " << usel.radius() << endl;

  vector<GCoord> bdd = usel.boundingBox();
  cout << "Bounding box is: ";
  cout << bdd[0] << " x " << bdd[1] << endl;

  GCoord c = p.centroid();
  cout << "The centroid for the PDB is at " << c << endl;

  c = usel.centroid();
  cout << "The centroid for the selection is at " << c << endl;

  AtomicGroup::Iterator iter(usel);
  int i;
  cout << "The first 5 atoms in the selection are...\n";
  for (i=0; i<5; i++)
    cout << *(iter()) << endl;

  PDB terminus = PDB::fromAtomicGroup(usel.subset(-1, 5));
  terminus.autoTerminate(false);
  cout << "\nThe last 5 are...\n";
  cout << terminus << endl;


  PDB split_ends = PDB::fromAtomicGroup(usel.subset(0, 5) + usel.subset(-1, 5));
  cout << "\nThe ends combined now...\n";
  cout << split_ends << endl;

  AtomicGroup residue = p.getResidue(usel[0]);
  residue.sort();
  cout << "\nThe first residue is:\n";
  cout << residue << endl;


}

