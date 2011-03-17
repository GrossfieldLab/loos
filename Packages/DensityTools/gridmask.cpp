/*
  gridmask.cpp


  (c) 2008 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry

  Given an int grid that represents picked blobs, use this as a mask
  against a double grid...
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
using namespace DensityTools;



int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "Usage- gridmask <edm_grid mask_grid >masked_edm_grid\n";
    exit(-1);
  }

  DensityGrid<int> mask;
  ifstream ifs(argv[1]);
  if (!ifs) {
    cerr << "Error - cannot open " << argv[1] << " for reading.\n";
    exit(-10);
  }

  ifs >> mask;

  DensityGrid<double> data;
  cin >> data;

  DensityGridpoint dims = data.gridDims();
  DensityGridpoint ddims = mask.gridDims();
  if (dims != ddims) {
    cerr << "Error - differing dimensions between mask and density grids.\n";
    exit(-10);
  }

  long k = dims[0] * dims[1] * dims[2];
  for (long i=0; i<k; i++)
    if (!mask(i))
      data(i) = 0.0;

  cout << data;
}
