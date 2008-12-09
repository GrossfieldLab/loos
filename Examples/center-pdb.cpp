/*
  center-pdb.cpp

  
  (c) 2008 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry

  Centers a PDB at the origin...  This is the LOOS equivalent of "hello world".
*/


#include <loos.hpp>

using namespace std;
using namespace loos;


int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "Usage- center-pdb pdb-file >output\n";
    exit(-1);
  }

  PDB pdb(argv[1]);
  pdb.centerAtOrigin();
  cout << pdb << endl;
}
