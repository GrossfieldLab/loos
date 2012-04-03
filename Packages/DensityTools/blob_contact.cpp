/*
  blob_contact.cpp

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

GCoord findMinCoord(const vector<GCoord>& coords) {
  GCoord min(numeric_limits<double>::max(),
             numeric_limits<double>::max(),
             numeric_limits<double>::max());

  for (vector<GCoord>::const_iterator i = coords.begin(); i != coords.end(); ++i)
    for (uint j=0; j<3; ++j)
      if ((*i)[j] < min[j])
        min[j] = (*i)[j];

  return(min);
}

GCoord findMaxCoord(const vector<GCoord>& coords) {
  GCoord max(numeric_limits<double>::min(),
             numeric_limits<double>::min(),
             numeric_limits<double>::min());

  for (vector<GCoord>::const_iterator i = coords.begin(); i != coords.end(); ++i)
    for (uint j=0; j<3; ++j)
      if ((*i)[j] > max[j])
        max[j] = (*i)[j];

  return(max);
}



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


vector<double> calculatePercentageContacts(const RealMatrix& M) {
  vector<long> sum(M.cols()-1, 0);

  for (uint j=0; j<M.rows(); ++j)
    for (uint i=1; i<M.cols(); ++i)
      sum[i-1] += M(j, i);

  vector<double> occ(M.cols()-1, 0.0);
  for (uint i=0; i<M.cols()-1; ++i)
    occ[i] = sum[i] / M.rows();

  return(occ);
}



int main(int argc, char *argv[]) {
  if (argc != 7) {
    cerr << "Usage- blob_contact model traj selection skip blobid distance <grid >out.asc 2>report.txt\n";
    cerr << "Note: requires an grid with blob ids (i.e. output from blobid)\n";
    exit(-1);
  }

  string hdr = invocationHeader(argc, argv);
  int k = 1;
  AtomicGroup model = createSystem(argv[k++]);
  pTraj traj = createTrajectory(argv[k++], model);
  AtomicGroup residue_subset = selectAtoms(model, argv[k++]);
  vector<AtomicGroup> residues = residue_subset.splitByResidue();
  uint skip = strtoul(argv[k++], 0, 10);
  int blobid = strtol(argv[k++], 0, 10);
  double distance = strtod(argv[k++], 0);

  DensityGrid<int> grid;
  cin >> grid;
  vector<GCoord> blob = findBlobCoords(grid, blobid);
  GCoord blobmin = findMinCoord(blob);
  GCoord blobmax = findMaxCoord(blob);

  
  ostringstream shdr;
  shdr << hdr << endl;
  shdr << boost::format("# Blob bounding box is %s x %s\n") % blobmin % blobmax;
  shdr << boost::format("# Blob has %d voxels\n") % blob.size();
  shdr << "# Residue list...\n";
  for (uint i=0; i<residues.size(); ++i)
    shdr << boost::format("# %d : %d %d %s %s\n")
      % i
      % residues[i][0]->id()
      % residues[i][0]->resid()
      % residues[i][0]->resname()
      % residues[i][0]->segid();
  
  uint n = traj->nframes();
  if (skip > 0) {
    n -= skip;
    traj->readFrame(skip-1);
  }

  RealMatrix M(n, residues.size()+1);
  uint t = 0;
    
  while (traj->readFrame()) {
    traj->updateGroupCoords(model);
    M(t, 0) = t + skip;
    vector<int> nearby = findResiduesNearBlob(blob, residues, distance);
    for (uint i=0; i<nearby.size(); ++i)
      M(t, i+1) = nearby[i];
    ++t;
  }
    

  writeAsciiMatrix(cout, M, shdr.str());


  cerr << "# " << hdr << endl;
  cerr << "# n\tresid\tatomid\tfractional contact\n";
  vector<double> fraction = calculatePercentageContacts(M);
  for (uint i=0; i<residues.size(); ++i) {
    cerr << boost::format("%d\t%d\t%d\t%f\n")
      % i
      % residues[i][0]->resid()
      % residues[i][0]->id()
      % fraction[i];
  }
}
