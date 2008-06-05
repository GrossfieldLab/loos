/*
  bounding.cpp

  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Displays the bounding box for a selection from a PDB...
*/



#include <loos.hpp>
#include <pdb.hpp>
#include <Selectors.hpp>
#include <Parser.hpp>

int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "Usage: " << argv[0] << " pdb-filename selection-string\n";
    exit(-1);
  }

  PDB pdb(argv[1]);
  Parser parsed(argv[2]);
  KernelSelector ksel(parsed.kernel());

  AtomicGroup subset = pdb.select(ksel);
  vector<GCoord> bdd = subset.boundingBox();
  cout << subset.size() << " atoms in subset.\n";
  cout << "Centroid at " << subset.centroid() << endl;
  cout << "Bounds: " << bdd[0] << " x " << bdd[1] << endl;

}
