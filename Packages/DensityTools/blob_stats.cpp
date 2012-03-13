/*
  blob_stats.cpp

  Gather statistics on blobs...
*/

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
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
#include <boost/format.hpp>
#include <boost/tuple/tuple.hpp>
#include <algorithm>
#include <sstream>
#include <limits>

#include <DensityGrid.hpp>

using namespace std;
using namespace loos;
using namespace loos::DensityTools;




int countBlobs(const DensityGrid<int>& grid) {
  DensityGridpoint dims = grid.gridDims();
  long k = dims[0] * dims[1] * dims[2];
  int n = 0;

  for (long i=0; i<k; i++)
    if (grid(i) > n)
      n = grid(i);

  return(n);
}


vector<int> sizeBlobs(const DensityGrid<int>& grid) {

  int n = countBlobs(grid);
  vector<int> sizes(n+1, 0);

  DensityGridpoint dims = grid.gridDims();
  long k = dims[0] * dims[1] * dims[2];
  
  for (long i=0; i<k; i++)
    sizes[grid(i)] += 1;
  
  return(sizes);
}


vector<GCoord> blobCentroids(const int n, const DensityGrid<int>& grid) {
  vector<GCoord> centers(n+1, GCoord(0,0,0));
  vector<int> counts(n+1, 0);

  DensityGridpoint dims = grid.gridDims();
  for (int k=0; k<dims[2]; k++)
    for (int j=0; j<dims[1]; j++)
      for (int i=0; i<dims[0]; i++) {
	DensityGridpoint point(i,j,k);
	GCoord u = grid.gridToWorld(point);
	int id = grid(point);

	centers[id] += u;
	counts[id] += 1;
      }

  for (int i=0; i<= n; i++)
    centers[i] /= counts[i];

  return(centers);
}


int main(int argc, char *argv[]) {
  DensityGrid<int> grid;

  if (argc != 1) {
    cerr <<
      "Usage- blob_stats <foo.grid\n"
      "\n"
      "Print out basic information about blobs in a grid (requires an integer grid)\n"
      "See also blobid, and pick_blob\n";
    exit(-1);
  }

  cin >> grid;
  cout << "Read in grid with dimensions " << grid.gridDims() << endl;
  cout << "Grid extents (real-space) is " << grid.minCoord() << " x " << grid.maxCoord() << endl;
  GCoord range = grid.maxCoord() - grid.minCoord();
  cout << "Grid range is " << range << endl;

  vector<int> sizes = sizeBlobs(grid);
  int n = sizes.size();
  vector<GCoord> centers = blobCentroids(n, grid);

  GCoord delta = grid.gridDelta();
  double voxel_volume = 1.0 / delta[0];
  for (int i=1; i<3; i++)
    voxel_volume *= (1.0 / delta[i]);

  cout << boost::format("Voxel volume = %8.6g\n") % voxel_volume;
  cout << boost::format("%6s %12s %12s\t%s\n")
    % "Id"
    % "Voxels"
    % "Size (in A^3)"
    % "Centroid (in A)";
  cout << boost::format("%-6s %-12s %-12s\t%s\n")
    % "------"
    % "------------"
    % "------------"
    % "------------------------------";



  for (uint i = 0; i < sizes.size(); ++i) {
    cout << boost::format("%6d %12d %12.6g\t") % i % sizes[i] % (sizes[i] * voxel_volume);
    cout << centers[i] << endl;
  }

}
