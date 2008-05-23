/*
  dcd-example2.cpp
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  An example of how to use the C++ DCD reader...  Doesn't really do anything
  useful...

*/

#include <ios>
#include <iostream>
#include <iomanip>

#include <loos.hpp>
#include <pdb.hpp>
#include <dcd.hpp>
#include <Selectors.hpp>


int main(int argc, char *argv[]) {

  // First, read in the PDB...
  PDB pdb(argv[1]);

  // Extract the non-solvent atoms...
  
  SolventSelector solvsel;
  NotSelector notsolvsel(solvsel);
  AtomicGroup nonsolv = pdb.select(notsolvsel);
  cout << "Found " << nonsolv.size() << " non-solvent atoms.\n";

  DCD dcd(argv[2]);

  int frameno = 0;
  while (dcd.readFrame()) {
    dcd.updateGroupCoords(nonsolv);
    cout << setw(6) << frameno++ << " = " << nonsolv.centroid() << endl;
  }

}
