/*
  grid2ascii

  (c) 2009 Tod D. Romo
*/





#include <loos.hpp>
#include <sgrid.hpp>

#include <boost/format.hpp>

using namespace std;
using namespace loos;
using namespace lab;


int main(int argc, char *argv[]) {
  SGrid<float> grid;

  cin >> grid;
  SGridpoint dim = grid.gridDims();
  GCoord min = grid.minCoord();
  GCoord max = grid.maxCoord();

  cout << boost::format("Read in grid of size %s\n") % dim;
  cout << boost::format("Grid range from %s x %s\n") % min % max;

  for (int k=0; k<dim.z(); ++k)
    for (int j=0; j<dim.y(); ++j)
      for (int i=0; i<dim.x(); ++i)
        cout << boost::format("(%d,%d,%d) = %f\n") % k % j % i % grid(k,j,i);
}
