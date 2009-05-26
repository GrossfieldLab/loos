/*
  bounding.cpp


  Displays the bounding box for a selection from a PDB...
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

using namespace std;
using namespace loos;

int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "Usage: " << argv[0] << " model-filename selection-string\n";
    exit(-1);
  }

  AtomicGroup model = createSystem(argv[1]);
  AtomicGroup subset = selectAtoms(model, argv[2]);
  vector<GCoord> bdd = subset.boundingBox();
  cout << subset.size() << " atoms in subset.\n";
  cout << "Centroid at " << subset.centroid() << endl;
  cout << "Bounds: " << bdd[0] << " x " << bdd[1] << endl;

  GCoord box = bdd[1];
  box.x() += fabs(bdd[0].x());
  box.y() += fabs(bdd[0].y());
  box.z() += fabs(bdd[0].z());
  cout << "Box: " << box << endl;
}
