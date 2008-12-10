/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Alan Grossfield
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



#include <iostream>
#include <string>
#include <vector>

#include <TimeSeries.hpp>

using namespace std;
using namespace loos;

int main(int argc, char *argv[]) {
    vector<double> tmp;
    for (int i=0; i<20; i++) {
        tmp.push_back(i);
    }

    TimeSeries foo;
    cout << foo.size() << endl;

    TimeSeries foo2(tmp);
    TimeSeries foo3(foo2);

    for (int i=0; i<20; i++) {
        cout << i << "  " << foo2[i] << "  " << foo3[i] << endl;
    }

    foo2 *= 5.0;

    TimeSeries foo4 = foo2 - foo3;
    TimeSeries foo5 = 1.0 + 2.0*foo3 - foo4;

    cout << endl;
    for (int i=0; i<20; i++) {
        cout << i << "  " << foo2[i] << "  " << foo4[i] << "  "
             << foo5[i] << endl;
    }

    cout << foo5.average() << "  "
         << foo5.variance() << "  "
         << foo5.stdev() << "  "
         << foo5.sterr() << "  "
         << endl;

    cout << "Test Block averaging" << endl;
    cout << 1 << "  " << foo5.block_var(1) << endl;
    cout << 2 << "  " << foo5.block_var(2) << endl;
    cout << 3 << "  " << foo5.block_var(3) << endl;
    cout << 4 << "  " << foo5.block_var(4) << endl;
    cout << 5 << "  " << foo5.block_var(5) << endl;
}
