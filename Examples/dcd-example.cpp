/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo
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
#include <iomanip>
#include <exception>
#include <vector>
#include <tr1/memory>

#include <stdio.h>
#include <assert.h>


#include <loos.hpp>
#include <utils.hpp>
#include <dcd.hpp>



int main(int argc, char *argv[]) {
  
  if (argc < 2 || argc > 3) {
    cerr << "Usage- " << argv[0] << " filename [flag]\n";
    exit(-1);
  }

  DCD dcd(argv[1]);

  cout << "DCD has " << dcd.natoms() << " in " << dcd.nsteps() << " steps.\n";
  cout << "Timestep is " << dcd.delta() << endl;
  if (dcd.hasCrystalParams())
    cout << "The DCD has crystal data.\n";

  if (argc == 2)
    exit(0);

  bool b;
  int n = 0;

  while ((b = dcd.readFrame())) {
    vector<dcd_double> xtal = dcd.crystalParams();
    
    vector<dcd_real> x = dcd.xcoords();
    vector<dcd_real> y = dcd.ycoords();
    vector<dcd_real> z = dcd.zcoords();

    int i;
    double c[3] = { 0.0, 0.0, 0.0 };
    float min[3] = { 1e30, 1e30, 1e30 };
    float max[3] = { -1e30, -1e30, -1e30 };

    for (i=0; i<dcd.natoms(); i++) {
      c[0] += x[i];
      c[1] += y[i];
      c[2] += z[i];
      
      if (x[i] < min[0])
	min[0] = x[i];
      if (y[i] < min[1])
	min[1] = y[i];
      if (z[i] < min[2])
	min[2] = z[i];

      if (x[i] > max[0])
	max[0] = x[i];
      if (y[i] > max[1])
	max[1] = y[i];
      if (z[i] > max[2])
	max[2] = z[i];

    }

    for (i=0; i<3; i++)
      c[i] /= dcd.natoms();

    // Remind me again why this is better than printf???
    cout << n << "\t" << setprecision(3) << setw(6) << xtal[0] << " x " << setw(6) << xtal[1] << " x " << setw(6) << xtal[2];
    cout << " = " << setprecision(3) << setw(8) << c[0] << " " << setw(8) << c[1] << " " << setw(8) << c[2];
    cout << "    -=> ";
    for (i=0; i<3; i++)
      cout << setw(8) << min[i];
    cout << " | ";
    for (i=0; i<3; i++)
      cout << setw(8) << max[i];
    cout << endl;

    ++n;
  }


}
