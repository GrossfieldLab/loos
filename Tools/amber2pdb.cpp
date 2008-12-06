#include <iostream>
#include <string>
#include <loos.hpp>

using namespace std;

int main(int argc, char *argv[]) {
  Amber file(argv[1], argv[2]);

  PDB pdb = PDB::fromAtomicGroup(file);
  cout << pdb;
}
