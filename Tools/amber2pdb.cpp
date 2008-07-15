#include <iostream>
#include <string>
#include <loos.hpp>


int main(int argc, char *argv[]) {
  Amber file(argv[1], argv[2]);

  PDB pdb = PDB::fromAtomicGroup(file);
  cout << pdb;
}
