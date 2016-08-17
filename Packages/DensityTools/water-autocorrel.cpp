/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2014, Tod D. Romo, Alan Grossfield
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


using namespace std;
using namespace loos;


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  if (argc == 1) {
    cerr << "Usage- " << argv[0] << " water_matrix [max-t] >output.asc\n";
    exit(-1);
  }

  int k = 1;
  string matname(argv[k++]);
  uint max_t = 0;
  if (k != argc)
    max_t = strtoul(argv[k++], 0, 10);
  
  Math::Matrix<int> M;
  cerr << "Reading matrix...\n";
  readAsciiMatrix(matname, M);
  uint m = M.rows();
  uint n = M.cols();

  if (max_t == 0)
    max_t = n/10;

  cerr << boost::format("Water matrix is %d x %d\n") % m % n;
  cerr << "Processing- ";
  vector< TimeSeries<double> > waters;
  for (uint j=0; j<m; ++j) {
    if (j % 250 == 0)
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

  uint nwaters = waters.size();
  cerr << boost::format(" done\nFound %d unique waters inside\n") % nwaters;
  cout << "# " << hdr << endl;
  for (uint j=0; j<max_t; ++j) {
    vector<double> tmp(nwaters, 0);
    for (uint i=0; i<nwaters; ++i)
      tmp[i] = waters[i][j];
    dTimeSeries dtmp(tmp);
    cout << j << '\t' << dtmp.average() << '\t' << dtmp.stdev() << '\t' << dtmp.sterr() << endl;
  }
}
