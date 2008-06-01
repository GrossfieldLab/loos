/*
  pdb2matlab.cpp
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Takes a PDB and a selection and an optional selection and writes out
  the coordinates to stdout in matlab format...

  The matrix is written in row-major order though, i.e.
  [ X_0 Y_0 Z_0 ;
    X_1 Y_1 Z_1 ;
       ...
    X_i Y_i Z_i ]

*/


#include <pdb.hpp>
#include <Parser.hpp>
#include <Selectors.hpp>


int main(int argc, char *argv[]) {
  
  if (argc  < 2 || argc > 3) {
    cerr << "Usage: " << argv[0] << " pdb-filename [selection string]" << endl;
    exit(-1);
  }

  PDB pdb(argv[1]);
  AtomicGroup atoms = pdb;

  if (argc > 2) {
    Parser parsed(argv[2]);
    KernelSelector parsed_selector(parsed.kernel());
    atoms = pdb.select(parsed_selector);
  }

  AtomicGroup::Iterator i(atoms);
  pAtom pa;

  cout << "A = [\n";
  while (pa = i())
    cout << pa->coords().x() << " " << pa->coords().y() << " " << pa->coords().z() << ";\n";
  cout << "];\n";

}
