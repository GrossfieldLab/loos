/*
  water-count.cpp

  (c) 2008 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry


   usage:
     water-count prefix >output.asc

   Description:
     Counts the number of waters that are inside the protein at each
     timepoint and writes this out as a vector.  In other words, it
     just sums the rows for each column of the water matrix.
*/



#include <loos.hpp>
#include <boost/format.hpp>


using namespace std;
using namespace loos;


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  if (argc != 2) {
    cerr << "Usage - water-count prefix >output.asc\n";
    exit(-1);
  }

  string prefix(argv[1]);

  Math::Matrix<double> V;
  readAsciiMatrix(prefix + ".vol", V);

  Math::Matrix<int> M;
  readAsciiMatrix(prefix + ".asc", M);
  uint m = M.rows();
  uint n = M.cols();

  if (V.rows() != M.cols()) {
    cerr << "ERROR - mismatch in volume and water matrix\n";
    exit(-10);
  }

  cout << "# " << hdr << endl;
  cout << "# t\tcount\tvolume\n";
  for (uint i=0; i<n; ++i) {
    long cnt = 0;
    for (uint j=0; j<m; ++j)
      cnt += M(j,i);
    cout << i << "\t" << cnt << "\t" << V(i,0) << endl;
  }
}
