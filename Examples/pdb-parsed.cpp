/*
  pdb-parsed.cpp
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School


  Basically same as pdb-example, but demonstrates a user-defined
  selector...
*/



#include <pdb.hpp>
#include <Parser.hpp>


int main(int argc, char *argv[]) {
  
  PDB p(argv[1]);

  cout << "Read in " << p.size() << " atoms from " << argv[1] << endl;

  Parser parsed(argv[2]);
  cout << "*** Virtual Machine Command STACK ***\n" << parsed.kernel() << endl;
  AtomicGroup  usel = p.select(parsed);
  
  cout << "There are " << usel.size() << " atoms in the selection.\n";
  cout << "The max radius is " << usel.radius() << endl;

  AtomicGroup::BoundingBox bdd = usel.boundingBox();
  cout << "Bounding box is:\n";
  cout << "\t" << bdd.min[0] << " <= x <= " << bdd.max[0] << endl;
  cout << "\t" << bdd.min[1] << " <= y <= " << bdd.max[1] << endl;
  cout << "\t" << bdd.min[2] << " <= z <= " << bdd.max[2] << endl;

  GCoord c = p.centroid();
  cout << "The centroid for the PDB is at " << c << endl;

  c = usel.centroid();
  cout << "The centroid for the selection is at " << c << endl;

  AtomicGroup::Iterator iter(usel);
  int i;
  cout << "The first 5 atoms in the selection are...\n";
  for (i=0; i<5; i++)
    cout << *(iter()) << endl;

  // Note the implicit conversion back to a PDB...
  PDB terminus = usel.subset(-1, 5);
  terminus.autoTerminate(false);
  cout << "\nThe last 5 are...\n";
  cout << terminus << endl;


  PDB split_ends = usel.subset(0, 5) + usel.subset(-1, 5);
  cout << "\nThe ends combined now...\n";
  cout << split_ends << endl;

  AtomicGroup residue = p.getResidue(usel[0]);
  residue.sort();
  cout << "\nThe first residue is:\n";
  cout << residue << endl;


}

