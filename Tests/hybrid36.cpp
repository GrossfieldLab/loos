
/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Alan Grossfield
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
#include <cmath>

using namespace std;
using namespace loos;


int ipow(int d, int n) {
  int p = 1;

  for (int i=0; i<n; ++i)
    p *= d;

  return(p);
}


int test(int n) {
  int max = 2*ipow(36,n);
  cerr << "* Testing n=" << n << " *\n";
  cerr << "Progress: ";
  for (int i=0; i<max; ++i) {

    if (i % 1000000 == 0)
      cerr << '.';

    try {
      string s = hybrid36AsString(i, n);
      int dd = parseStringAsHybrid36(s);
      if (dd != i) {
        cout << "***ERROR***\n";
        cout << "i = " << i << endl;
        cout << "s = ]" << s << "[\n";
        cout << "dd = " << dd << endl;
        exit(-1);
      }
    }
    catch (...) {
      cerr << "done\n";
      return(i);
    }
  }

  return(-1);
}


int main(int argc, char *argv[]) {

  int n = test(4);
  if (n != 2436112) {
    cerr << "FAILED: n=" << n << endl;
    exit(-1);
  }

  n = test(5);
  if (n != 87440032) {
    cerr << "Failed: n=" << n << endl;
    exit(-1);
  }
}
