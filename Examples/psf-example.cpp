/*
  psf-example.cpp
  (c) 2008 Alan Grossfield

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School


  Some examples of using the PDB/AtomicGroup classes...
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
    AtomicGroup::BoundingBox bdd = water.boundingBox();
    cout << "Bounding box for the water is:\n";
    cout << "\t" << bdd.min[0] << " <= x <= " << bdd.max[0] << endl;
    cout << "\t" << bdd.min[1] << " <= y <= " << bdd.max[1] << endl;
    cout << "\t" << bdd.min[2] << " <= z <= " << bdd.max[2] << endl;
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

