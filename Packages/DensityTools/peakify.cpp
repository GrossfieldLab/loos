/*
  Peakify.cpp

  (c) 2009 Tod D. Romo, Grossfield Lab, URMC


  Creates a PDB representing peaks in the grid...

  Usage:  peakify threshold <input.grid >output.pdb
*/



#include <loos.hpp>
#include <sgrid.hpp>
#include <sgrid_utils.hpp>

using namespace std;
using namespace loos;
using namespace lab;


int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "Usage- peakify threshold <grid >pdb\n";
    exit(-1);
  }

  double thresh = strtod(argv[1], 0);

  string hdr = invocationHeader(argc, argv);
  SGrid<double> grid;
  cin >> grid;

  cerr << "Read in grid " << grid.gridDims() << "\n";

  SGrid<int> blobs(grid.minCoord(), grid.maxCoord(), grid.gridDims());

  PDB pdb;
  vector<GCoord> peaks = findPeaks(grid, Threshold<double>(thresh));
  for (int i = 0; i < peaks.size(); ++i) {
    pAtom atom(new Atom(i+1, "UNK", peaks[i]));
    atom->resid(i+1);
    atom->resname("UNK");
    atom->segid("BLOB");

    pdb.append(atom);
  }

  pdb.remarks().add(hdr);
  cout << pdb;
}
