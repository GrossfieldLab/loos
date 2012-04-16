/*
  Assign frames of a trajectory to bins (determined by reference structures)
*/



/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2011, Tod D. Romo
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

#include "fid-lib.hpp"

using namespace std;
using namespace loos;


string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\tConstruct a structural histogram given a set of fiducial structures\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tThis tool assigns the frames of a trajectory to the closest bin based\n"
    "on the fiducial structures given.  See Lyman & Zuckerman,\n"
    "J Phys Chem B (2007) 111:12876-12882 for more details.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\tassign_frames model.pdb simulation.dcd all 'name == \"CA\"' zuckerman.dcd >assignments.asc\n"
    "This example assigns all frames in simulation.dcd using the fiducials stored in zuckerman.dcd,\n"
    "writing the assignments to assignments.asc.\n"
    "\n"
    "NOTES\n"
    "\tThe selection used here must match that given to ufidpick\n"
    "SEE ALSO\n"
    "\tufidpick, effsize.pl, hierarchy\n";

  return(msg);
}


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  if (argc != 6) {
    cerr << "Usage - " << argv[0] << " model trajectory range selection fiducials.dcd >assignments.asc\n";
    fullHelpMessage();
    exit(-1);
  }
  
  int k = 1;
  AtomicGroup model = createSystem(argv[k++]);
  pTraj traj = createTrajectory(argv[k++], model);
  string range(argv[k++]);
  string selection(argv[k++]);
  
  AtomicGroup subset = selectAtoms(model, selection);
  AtomicGroup ref_model = subset.copy();
  ref_model.renumber();
  pTraj fiducials = createTrajectory(argv[k++], ref_model);

  vecUint frames;
  if (range == "all")
    for (uint i=0; i<traj->nframes(); ++i)
      frames.push_back(i);
  else
    frames = parseRangeList<uint>(range);

  vecGroup refs;
  cerr << "Reading fiducials...\n";
  readTrajectory(refs, ref_model, fiducials);
  cerr << "Read in " << refs.size() << " fiducials.\n";
  cerr << "Assigning...\n";
  vecUint assigned = assignStructures(subset, traj, frames, refs);
  cout << "# " << hdr << endl;
  copy(assigned.begin(), assigned.end(), ostream_iterator<uint>(cout, "\n"));

}
