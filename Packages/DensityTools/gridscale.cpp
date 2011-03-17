/*
  gridscale.cpp
  (c) 2009 Tod D. Romo, Grossfield Lab, URMC

  Applies a constant scaling to a grid...
*/



#include <loos.hpp>
#include <sgrid.hpp>

using namespace std;
using namespace loos;
using namespace lab;

int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  if (argc != 2) {
    cerr << "Usage- gridscale scale-value <in-grid >out-grid\n";
    exit(0);
  }

  double konst = strtod(argv[1], 0);

  SGrid<double> grid;
  cin >> grid;
  grid.scale(konst);
  grid.addMetadata(hdr);
  cout << grid;
}
  
