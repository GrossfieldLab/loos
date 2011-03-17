/*
  gridslice.cpp


  (c) 2008 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry

  Takes a double grid and extracts a plane from it as a matrix...
*/


#include <loos.hpp>
#include <boost/format.hpp>
#include <boost/tuple/tuple.hpp>
#include <sgrid.hpp>

using namespace std;
using namespace loos;
using namespace lab;

typedef Math::Matrix<double, Math::RowMajor> Matrix;

void invalidIndex(int i) {
  cerr << "ERROR - invalid plane index " << i << endl;
  exit(-1);
}



int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "Usage - gridslice [i|j|k] index <grid >matrix\n";
    exit(-1);
  }

  string hdr = invocationHeader(argc, argv);

  string plane(argv[1]);
  int idx = atoi(argv[2]);

  SGrid<double> grid;
  cin >> grid;
  SGridpoint dims = grid.gridDims();
  cerr << boost::format("Grid dimensions are %d x %d x %d (i x j x k)\n") % dims[0] % dims[1] % dims[2];
  if (plane == "k") {

    if (idx > dims[2])
      invalidIndex(idx);
    Matrix M(dims[0]+1, dims[1]+1);
    for (int j=0; j<dims[1]; ++j)
      for (int i=0; i<dims[0]; ++i)
        M(j,i) = grid(idx,j,i);

    writeAsciiMatrix(cout, M, hdr);

  } else if (plane == "j") {

    if (idx > dims[1])
      invalidIndex(idx);
    Matrix M(dims[0]+1, dims[2]+1);
    for (int k=0; k<dims[2]; ++k)
      for (int i=0; i<dims[0]; ++i)
        M(k,i) = grid(k,idx,i);
    
    writeAsciiMatrix(cout, M, hdr);

  } else if (plane == "i") {

    if (idx > dims[0])
      invalidIndex(idx);
    Matrix M(dims[1]+1, dims[2]+1);
    for (int k=0; k<dims[2]; ++k)
      for (int j=0; j<dims[1]; ++j)
        M(k,j) = grid(k,j,idx);

    writeAsciiMatrix(cout, M, hdr);

  } else {
    cerr << "ERROR - unknown plane '" << plane << "'\n";
    exit(-1);
  }

}
