// (c) 2014 Tod D. Romo, Grossfield Lab, URMC
//

#include <loos.hpp>


using namespace std;
using namespace loos;


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  int k = 1;
  string matname(argv[k++]);
  uint max_t = 0;
  if (k != argc)
    max_t = strtoul(argv[k++], 0, 10);
  
  Math::Matrix<int> M;
  readAsciiMatrix(matname, M);
  uint m = M.rows();
  uint n = M.cols();

  if (max_t == 0)
    max_t = n/2;

  cerr << boost::format("Water matrix is %d x %d\n") % m % n;
  cerr << "Processing- ";
  vector< TimeSeries<double> > waters;
  for (uint j=0; j<m; ++j) {
    if (j % 50 == 0)
      cerr << '.';

    vector<double> tmp(n, 0);
    bool flag = false;
    for (uint i=0; i<n; ++i) {
      tmp[i] = M(j, i);
      if (tmp[i])
	flag = true;
    }
    if (flag) {
      TimeSeries<double> wtmp(tmp);
      waters.push_back(wtmp.correl(max_t));
    }
  }
  cerr << boost::format(" done\nFound %d unique waters inside\n") % waters.size();
  cout << "# " << hdr << endl;
  for (uint j=0; j<max_t; ++j) {
    vector<double> tmp(m, 0);
    for (uint i=0; i<m; ++i)
      tmp[i] = waters[i][j];
    dTimeSeries dtmp(tmp);
    cout << j << '\t' << dtmp.average() << '\t' << dtmp.stdev() << '\t' << dtmp.sterr();
    for (uint i=0; i<m; ++i)
      cout << '\t' << waters[i][j];
    cout << endl;
  }
}
