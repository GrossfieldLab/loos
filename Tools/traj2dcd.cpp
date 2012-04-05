/*
  traj2dcd


  Converts a LOOS-supported format to a DCD

  Usage:

    traj2dcd model-file trajectory-file dcd-name

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
#include <boost/format.hpp>

using namespace std;
using namespace loos;



string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\tConvert trajectory into DCD format\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tConvert any LOOS-supported trajectory format into a DCD trajectory.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\ttraj2dcd model.gro simulation.xtc simulation.dcd\n"
    "Convert the GROMACS XTC trajectory into the DCD format.\n"
    "\n"
    "SEE ALSO\n"
    "\tsubsetter, merge-traj, recenter-traj, reimage-by-molecule\n";

  return(msg);
}



int main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << "Usage - traj2dcd model trajectory dcd\n";
    cerr << fullHelpMessage();
    exit(-1);
  }

  AtomicGroup model = createSystem(argv[1]);
  pTraj traj = createTrajectory(argv[2], model);
  uint n = traj->nframes();

  DCDWriter dcd(argv[3]);
  dcd.setHeader(model.size(), n, 1e-3, traj->hasPeriodicBox());
  dcd.setTitle(invocationHeader(argc, argv));
  dcd.writeHeader();

  cerr << boost::format("There are %d atoms and %d frames.\n") % model.size() % n;

  cerr << "Processing - ";
  cerr.flush();

  for (uint i=0; i<n; ++i) {
    if (i % 250 == 0)
      cerr << '.';
    traj->readFrame(i);
    traj->updateGroupCoords(model);
    dcd.writeFrame(model);
  }
  
  cerr << " done\n";
}

