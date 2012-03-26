/*
  near_blob.cpp

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



vector<GCoord> findBlobCoords(DensityGrid<int>& grid, const int blobid) {
  DensityGridpoint griddims = grid.gridDims();
  
  vector<GCoord> coords;
  for (int k=0; k<griddims.z(); ++k)
    for (int j=0; j<griddims.y(); ++j)
      for (int i=0; i<griddims.x(); ++i) {
        DensityGridpoint grid_coord(i, j, k);
        if (grid(grid_coord) == blobid)
          coords.push_back(grid.gridToWorld(grid_coord));
      }

  return(coords);
}


vector<int> findResiduesNearBlob(const vector<GCoord>& blob, const vector<AtomicGroup>& residues, const double threshold) {
  
  double thresh2 = threshold * threshold;
  vector<int> nearby(residues.size(), 0);

  for (uint k=0; k<residues.size(); ++k) {
    bool flag = false;
    for (uint j=0; j<residues[k].size() && !flag; ++j) {
      GCoord c = residues[k][j]->coords();

      for (uint i=0; i<blob.size(); ++i) {
        double d = c.distance2(blob[i]);
        if (d <= thresh2) {
          flag = true;
          break;
        }
      }
    }
    nearby[k] = flag;
  }

  return(nearby);
}



int main(int argc, char *argv[]) {
  if (argc != 6) {
    cerr << "Usage- near_blob model traj selection blobid distance <grid >out.asc\n";
    exit(-1);
  }

  string hdr = invocationHeader(argc, argv);
  int k = 1;
  AtomicGroup model = createSystem(argv[k++]);
  pTraj traj = createTrajectory(argv[k++], model);
  AtomicGroup residue_subset = selectAtoms(model, argv[k++]);
  vector<AtomicGroup> residues = residue_subset.splitByResidue();
  int blobid = strtol(argv[k++], 0, 10);
  double distance = strtod(argv[k++], 0);

  DensityGrid<int> grid;
  cin >> grid;
  vector<GCoord> blob = findBlobCoords(grid, blobid);
  
  ostringstream shdr;
  shdr << hdr << endl;
  shdr << "# Residue list...\n";
  for (uint i=0; i<residues.size(); ++i)
    shdr << boost::format("# %d : %d %d %s %s\n")
      % i
      % residues[i][0]->id()
      % residues[i][0]->resid()
      % residues[i][0]->resname()
      % residues[i][0]->segid();
  
  RealMatrix M(traj->nframes(), residues.size()+1);
  uint t = 0;
    
  while (traj->readFrame()) {
    M(t, 0) = t;
    vector<int> nearby = findResiduesNearBlob(blob, residues, distance);
    for (uint i=0; i<nearby.size(); ++i)
      M(t, i+1) = nearby[i];
    ++t;
  }
    

  writeAsciiMatrix(cout, M, shdr.str());
}
