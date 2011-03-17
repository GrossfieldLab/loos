/*
  blob_stats.cpp


  (c) 2008 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry

  Gather statistics on blobs...
*/


#include <loos.hpp>
#include <boost/format.hpp>
#include <boost/tuple/tuple.hpp>
#include <algorithm>
#include <sstream>
#include <limits>

#include <sgrid.hpp>

using namespace std;
using namespace loos;




int countBlobs(const lab::SGrid<int>& grid) {
  lab::SGridpoint dims = grid.gridDims();
  long k = dims[0] * dims[1] * dims[2];
  int n = 0;

  for (long i=0; i<k; i++)
    if (grid(i) > n)
      n = grid(i);

  return(n);
}


vector<int> sizeBlobs(const lab::SGrid<int>& grid) {

  int n = countBlobs(grid);
  vector<int> sizes(n+1, 0);

  lab::SGridpoint dims = grid.gridDims();
  long k = dims[0] * dims[1] * dims[2];
  
  for (long i=0; i<k; i++)
    sizes[grid(i)] += 1;
  
  return(sizes);
}


vector<GCoord> blobCentroids(const int n, const lab::SGrid<int>& grid) {
  vector<GCoord> centers(n+1, GCoord(0,0,0));
  vector<int> counts(n+1, 0);

  lab::SGridpoint dims = grid.gridDims();
  for (int k=0; k<dims[2]; k++)
    for (int j=0; j<dims[1]; j++)
      for (int i=0; i<dims[0]; i++) {
	lab::SGridpoint point(i,j,k);
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
  lab::SGrid<int> grid;

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



  for (int i = 0; i < sizes.size(); ++i) {
    cout << boost::format("%6d %12d %12.6g\t") % i % sizes[i] % (sizes[i] * voxel_volume);
    cout << centers[i] << endl;
  }

}
