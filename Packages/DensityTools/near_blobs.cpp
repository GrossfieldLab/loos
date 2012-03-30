/*
  near_blobs.cpp

  Find residues within a given distance of a blob

*/

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2012, Tod D. Romo, Alan Grossfield
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
#include <DensityTools.hpp>

using namespace std;
using namespace loos;
using namespace loos::DensityTools;


typedef vector<GCoord>        vCoords;
typedef vector<vCoords>      vvCoords;
typedef vector<AtomicGroup>    vGroup;



vvCoords separateBlobs(const DensityGrid<int>& grid) {

  int max_blobid = 0;

  for (long i = 0; i < grid.size(); ++i)
    if (grid(i) > max_blobid)
      max_blobid = grid(i);

  vvCoords blobs(max_blobid);

  DensityGridpoint dims = grid.gridDims();
  for (int k=0; k<dims.z(); ++k)
    for (int j=0; j<dims.y(); ++j)
      for (int i=0; i<dims.x(); ++i) {
        DensityGridpoint p(i, j, k);
        if (grid(p)) {
          int id = grid(p);
          GCoord c = grid.gridToWorld(p);
          blobs[id-1].push_back(c);
        }
      }
  
  return(blobs);
}



vector<uint> findBlobsNearResidue(const vvCoords& blobs, const AtomicGroup& residue, const double dist) {
  vector<uint> blobids;

  double d2 = dist * dist;
  for (uint k=0; k<blobs.size(); ++k) {
    bool flag = true;

    for (uint j=0; j<residue.size() && flag; ++j) {
      GCoord c = residue[j]->coords();

      for (uint i=0; i<blobs[k].size() && flag; ++i)
        if (c.distance2(blobs[k][i]) <= d2)
          flag = false;

    }

    if (!flag)
      blobids.push_back(k);

  }

  return(blobids);
}




int main(int argc, char *argv[]) {
  
  if (argc != 4) {
    cerr << "Usage- near_blobs model selection distance <blobs.grid\n";
    cerr << "NOTE: grid must have IDs (i.e. output from blobid)\n";
    exit(-1);
  }

  string hdr = invocationHeader(argc, argv);

  int k = 1;
  AtomicGroup model = createSystem(argv[k++]);
  AtomicGroup subset = selectAtoms(model, argv[k++]);
  double distance = strtod(argv[k++], 0);

  DensityGrid<int> the_grid;
  cin >> the_grid;

  vvCoords blobs = separateBlobs(the_grid);
  vGroup residues = subset.splitByResidue();

  cout << "# " << hdr << endl;
  cout << "# Atomid Resid Resname Segid Bloblist...\n";
  for (uint i=0; i<residues.size(); ++i) {
    vector<uint> ids = findBlobsNearResidue(blobs, residues[i], distance);
    if (ids.size() == 0)
      continue;
    cout << boost::format("%d\t%d\t%s\t%s\t")
      % residues[i][0]->id()
      % residues[i][0]->resid()
      % residues[i][0]->resname()
      % residues[i][0]->segid();
    for (uint j=0; j<ids.size(); ++j)
      cout << ids[j]+1 << (j == ids.size()-1 ? "" : ",");
    cout << endl;
  }
}
