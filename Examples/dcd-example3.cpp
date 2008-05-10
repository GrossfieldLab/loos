/*
  dcd-example3.cpp
  (c) 2008 Tod D. Romo and Alan Grossfield

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  An example of how to use the C++ DCD reader...  Doesn't really do anything
  useful...

  Same as dcd-example2, except this time we read a psf to define the system
  instead of a pdb

*/

#include <ios>
#include <iostream>
#include <iomanip>

#include <loos.hpp>
#include <psf.hpp>
#include <dcd.hpp>


// Select for non-solvent atoms...
struct NotSolvSelector : public AtomSelector {
  bool operator()(const pAtom& atom) const {
    return(!(atom->segid() == "SOLV" || atom->segid() == "BULK"));
  }
};



int main(int argc, char *argv[]) {

  // First, read in the PSF...
  PSF psf(argv[1]);

  // Extract the non-solvent atoms...
  NotSolvSelector ns;
  AtomicGroup nonsolv = psf.select(ns);
  cout << "Found " << nonsolv.size() << " non-solvent atoms.\n";

  DCD dcd(argv[2]);

  int frameno = 0;
  while (dcd.readFrame()) {
    dcd.updateGroupCoords(nonsolv);
    cout << setw(6) << frameno++ << " = " << nonsolv.centroid() << endl;
  }

}
