/*
  Peakify.cpp

  Creates a PDB representing peaks in the grid...

  Usage:  peakify threshold <input.grid >output.pdb
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo, Alan Grossfield
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester

  This package (LOOS) is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation under version 3 of the License.

  This package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <loos.hpp>
#include <DensityGrid.hpp>
#include <GridUtils.hpp>

using namespace std;
using namespace loos;
using namespace loos::DensityTools;


int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "Usage- peakify threshold <grid >pdb\n";
    exit(-1);
  }

  double thresh = strtod(argv[1], 0);

  string hdr = invocationHeader(argc, argv);
  DensityGrid<double> grid;
  cin >> grid;

  cerr << "Read in grid " << grid.gridDims() << "\n";

  DensityGrid<int> blobs(grid.minCoord(), grid.maxCoord(), grid.gridDims());

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
