/*
  pdb-example.cpp
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School


  Some examples of using the PDB/AtomicGroup classes...
*/


#include <pdb.hpp>

int main(int argc, char *argv[]) {

  PDB pdb(argv[1]);
  pdb.centerAtOrigin();

  cout << pdb << endl;
}

