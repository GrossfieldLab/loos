/*
  gridinfo.cpp
  (c) 2009 Tod D. Romo, Grossfield Lab, URMC


  Just dump the grid header info...
*/



#include <loos.hpp>
#include <sgrid.hpp>


using namespace std;
using namespace lab;



int main(int argc, char *argv[]) {

  SGrid<double> grid;
  if (argc == 2) {
    ifstream ifs(argv[1]);
    ifs >> grid;
  } else
    cin >> grid;

  loos::GCoord min = grid.minCoord();
  loos::GCoord max = grid.maxCoord();
  SGridpoint dim = grid.gridDims();

  cout << "Grid = " << min << " x " << max << " @ " << dim << endl;
  cout << "Resolution = " << (max.x() - min.x()) / dim.x() << endl;
  cout << "Metadata: ";

  SMetaData meta = grid.metadata();
  if (meta.empty())
    cout << "none\n";
  else {
    cout << endl;
    copy(meta.begin(), meta.end(), ostream_iterator<SMetaData::value_type>(cout, "\n"));
  }
}
