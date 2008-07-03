#include <loos.hpp>

int main(int argc, char *argv[]) {
  
  if (argc != 2) {
    cerr << "Usage- pdb-echo <pdb-file>\n";
    exit(-1);
  }

  PDB pdb(argv[1]);
  cout << pdb;
}
