/*
  Gridgauss

  (c) 2009 Tod D. Romo
*/


#include <loos.hpp>
#include <boost/format.hpp>
#include <sgrid.hpp>
#include <sgrid_utils.hpp>

#include "banal-lib.hpp"

using namespace std;
using namespace loos;
using namespace lab;


int main(int argc, char *argv[]) {

  if (argc != 3) {
    cerr << "Usage- gridgauss width sigma <grid >output\n";
    exit(-1);
  }
  
  string hdr = invocationHeader(argc, argv);
  int width = atoi(argv[1]) - 1;
  double sigma = strtod(argv[2], 0);

  vector<double> kernel = banal::gaussian1d(width, sigma);
  cerr << "Kernel (" << kernel.size() << "): ";
  copy(kernel.begin(), kernel.end(), ostream_iterator<double>(cerr, ","));
  cerr << endl;

  SGrid<double> grid;
  cin >> grid;
  gridConvolve(grid, kernel);

  grid.addMetadata(hdr);
  cout << grid;
}


  
