/*
  pdb-example.cpp
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School


  Some examples of using the PDB/AtomicGroup classes...
*/


#include <utils.hpp>
#include <pdb.hpp>
#include <Selectors.hpp>


int main(int argc, char *argv[]) {
  
  PDB p(argv[1]);

  cout << "Read in " << p.size() << " atoms from " << argv[1] << endl;
  if (remarksHasBox(p.remarks())) {
    GCoord box = boxFromRemarks(p.remarks());
    cout << "Alan-box detected: " << box << endl;
  } else
    cout << "No Alan-box detected!\n";

  CAlphaSelector casel;
  AtomicGroup cas = p.select(casel);
  
  cout << "There are " << cas.size() << " CAs.\n";
  cout << "The max radius for CAs is " << cas.radius() << endl;

  SolventSelector wasel;
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

  NotSelector notwatsel(wasel);
  AtomicGroup notwatgrp = p.select(notwatsel);
  cout << "There are " << notwatgrp.numberOfResidues() << " residues that are non-solvent.\n";

  GCoord c = p.centroid();
  cout << "The centroid for the PDB is at " << c << endl;

  AtomicGroup::Iterator iter(cas);
  int i;
  cout << "The first 5 CAs are...\n";
  for (i=0; i<5; i++)
    cout << *(iter()) << endl;

  PDB terminus = PDB::fromAtomicGroup(cas.subset(-1, 5));
  terminus.autoTerminate(false);
  cout << "\nThe last 5 CA's are...\n";
  cout << terminus << endl;


  PDB split_ends = PDB::fromAtomicGroup(cas.subset(0, 5) + cas.subset(-1, 5));
  cout << "\nThe ends combined now...\n";
  cout << split_ends << endl;

  AtomicGroup residue = p.getResidue(cas[0]);
  residue.sort();
  cout << "\nThe first residue is:\n";
  cout << residue << endl;

  HeavyAtomSelector heav_sel;
  HydrogenSelector hyd_sel;
  AtomicGroup hydrogens = residue.select(hyd_sel);
  AtomicGroup heavy_atoms = residue.select(heav_sel);

  cout << "Hydrogens   " << hydrogens.size() << endl;
  cout << hydrogens << endl;
  cout << "Heavy   " << heavy_atoms.size() << endl;
  cout << heavy_atoms << endl;

}

