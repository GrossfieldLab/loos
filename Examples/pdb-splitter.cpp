/*
  pdb-parsed.cpp
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School


*/



#include <pdb.hpp>
#include <Parser.hpp>
#include <Selectors.hpp>


int main(int argc, char *argv[]) {
  
  PDB p(argv[1]);

  cout << "Read in " << p.size() << " atoms from " << argv[1] << endl;

  Parser parsed(argv[2]);
  KernelSelector parsed_sel(parsed.kernel());
  AtomicGroup  usel = p.select(parsed_sel);
  
  cout << "There are " << usel.size() << " atoms in the selection.\n";
  vector<AtomicGroup> grps = usel.splitByUniqueSegid();

  cout << "There are " << grps.size() << " groups selected.\n";

  cout << grps[0] << endl;

}

