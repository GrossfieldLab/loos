/*
  traj2matlab.cpp

  Matrix written is in row-major order
  [ x x x ... ]
  [ y y y ... ]
  [ z z z ... ]
  [ . . . ... ]
  [ . . . ... ]
  [ . . . ... ]
*/




/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo
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

typedef Math::Matrix<double, Math::RowMajor>   Matrix;


string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\tConvert a trajectory into an ASCII matrix representation\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tThis tool will extract a subset of atoms from a trajectory\n"
    "and write it out as an ASCII matrix suitable for reading into\n"
    "octave and matlab.  Each frame of the trajectory becomes a column\n"
    "in the 3NxT matrix where T is the number of frames and N is the number\n"
    "of atoms.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\ttraj2matlab model.pdb simulation.dcd 'segid == \"PROT\" && !hydrogen' >M.asc\n"
    "This writes out all non-hydrogen atoms in the PROT segment to M.asc.\n"
    "\n"
    "NOTES\n"
    "\tA PDB is both a model and a single-frame trajectory.  A single model can therefore\n"
    "be converted by using the same file for both the model and the trajectory, i.e.\n"
    "\t\ttraj2matlab model.pdb model.pdb 'all' >model.asc\n"
    "SEE ALSO\n"
    "\tsvd\n";

  return(msg);
}



int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  if (argc  != 4) {
    cerr << "Usage: " << argv[0] << " model trajectory selection" << endl;
    cerr << fullHelpMessage();
    exit(-1);
  }

  AtomicGroup model = createSystem(argv[1]);
  pTraj traj = createTrajectory(argv[2], model);
  AtomicGroup subset = selectAtoms(model, argv[3]);

  Matrix A(subset.size() * 3, traj->nframes());
  for (uint i = 0; i<traj->nframes(); ++i) {
    traj->readFrame(i);
    traj->updateGroupCoords(subset);
    
    for (uint j=0, k=0; j<subset.size(); ++j) {
      GCoord c = subset[j]->coords();
      A(k++,i) = c.x();
      A(k++,i) = c.y();
      A(k++,i) = c.z();
    }
  }

  writeAsciiMatrix(cout, A, hdr);
}
