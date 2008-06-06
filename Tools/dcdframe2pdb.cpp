/*
  dcdframe2pdb

  dcdframe2pdb pdb dcd frameno >output
*/


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





#include <loos.hpp>
#include <dcd.hpp>
#include <pdb.hpp>


int main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << argv[0] << " pdbfile dcdfile frameno\n";
    exit(-1);
  }

  PDB pdb(argv[1]);
  DCD dcd(argv[2]);
  int frame = atoi(argv[3]);


  bool b = dcd.readFrame(frame);
  if (!b) {
    cerr << "Could not read frame " << frame << " from DCD " << argv[2] << endl;
    exit(-2);
  }


  dcd.updateGroupCoords(pdb);

  cout << pdb;

}


