/*
  dcdframe2pdb
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School


  dcdframe2pdb pdb dcd frameno >output
*/


#include <loos.hpp>
#include <dcd.hpp>
#include <pdb.hpp>


int main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << argv[0] << " pdbfile dcdfile frameno\n";
    exit(-1);
  }

  PDB pdb(argv[1]);
  DCD dcd(argv[2]);
  int frame = atoi(argv[3]);


  bool b = dcd.readFrame(frame);
  if (!b) {
    cerr << "Could not read frame " << frame << " from DCD " << argv[2] << endl;
    exit(-2);
  }


  dcd.updateGroupCoords(pdb);

  cout << pdb;

}


