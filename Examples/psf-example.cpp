/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
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




#include <psf.hpp>


struct CASelector : public AtomSelector {
  bool operator()(const pAtom& atom) const {
    return(atom->name() == "CA");
  }
};



struct SolvSelector : public AtomSelector {
  bool operator()(const pAtom& atom) const  {
    return(atom->segid() == "SOLV" || atom->segid() == "BULK");
  }
};


int main(int argc, char *argv[]) {
  
  PSF p(argv[1]);

  cout << "Read in " << p.size() << " atoms from " << argv[1] << endl;

  CASelector casel;
  AtomicGroup cas = p.select(casel);
  
  cout << "There are " << cas.size() << " CAs.\n";
  cout << "The max radius for CAs is " << cas.radius() << endl;

  SolvSelector wasel;
  AtomicGroup water = p.select(wasel);

  int nwater = water.numberOfResidues();
  cout << "There are " << nwater << " waters.\n";
  if (nwater > 0) {
    vector<GCoord> bdd = water.boundingBox();
    cout << "Bounding box for the water is: ";
    cout << bdd[0] << " x " << bdd[1] << endl;
  }

  GCoord c = p.centroid();
  cout << "The centroid for the PDB is at " << c << endl;

  AtomicGroup::Iterator iter(cas);
  int i;
  cout << "The first 5 CAs are...\n";
  for (i=0; i<5; i++)
    cout << *(iter()) << endl;

  AtomicGroup residue = p.getResidue(cas[0]);
  residue.sort();
  cout << "\nThe first residue is:\n";
  cout << residue << endl;

  cout << "Test groupFromID";
  pAtom pa = residue.getAtom(0);
  cout << "Atom: " << *pa << endl;
  vector<int> bondIDs = pa->getBonds();
  for (unsigned int i=0; i<bondIDs.size(); i++) { cout << bondIDs[i] << "  ";}
  cout << endl;
  AtomicGroup bonded = p.groupFromID(bondIDs);
  cout << bonded << endl;

}

