/*
  ufidpick
  
  Picks uniform fiducial structures for a structural histogram, a la Lyman &
  Zuckerman, J Phys Chem B (2007) 111:12876-12882


  Usage- ufidpick model trajectory selection output-name probability [seed]

*/



/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2010, Tod D. Romo
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
    "\tPick fiducial structures for a structural histogram\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tThis tool picks structures from a trajectory for use as fiducials in a\n"
    "structural histogram.  They are picked using bins with a uniform probability.  For\n"
    "more details, see Lyman & Zuckerman, J Phys Chem B (2007) 111:12876-12882.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\tufidpick model.pdb simulation.dcd all 'name == \"CA\"' zuckerman 0.1 >ufidpick.log\n"
    "This example uses bins with a probability of 0.1 (i.e. 10 bins), using only\n"
    "the alpha-carbons.  The output files include a log of what structures were \n"
    "picked, stored in ufidpick.log, as well as a trajectory containing just the\n"
    "fiducial structures in zuckerman.dcd and the corresponding model file in zuckerman.pdb\n"
    "\n"
    "SEE ALSO\n"
    "\tassign_frames, hierarchy, effsize.pl, neff\n";

  return(msg);
}


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  if (argc < 7 || argc > 8) {
    cerr << "Usage - " << argv[0] << " model trajectory range|all selection output-name cutoff [seed]\n";
    cerr << fullHelpMessage();
    exit(-1);
  }

  int opti = 1;
  AtomicGroup model = createSystem(argv[opti++]);
  model.clearBonds();

  pTraj traj = createTrajectory(argv[opti++], model);
  string range(argv[opti++]);
  string selection(argv[opti++]);
  AtomicGroup subset = selectAtoms(model, selection);
  string outname(argv[opti++]);
  double cutoff = strtod(argv[opti++], 0);

  uint seed;
  if (opti < argc) {
    seed = strtoul(argv[opti++], 0, 10);
    rng_singleton().seed(static_cast<unsigned int>(seed));  // This addresses a bug in earlier BOOSTs
  } else
    seed = randomSeedRNG();


  cout << "# " << hdr << endl;
  cout << "# seed = " << seed << endl;


  vector<uint> source_frames;
  if (range == "all")
    for (uint i=0; i<traj->nframes(); ++i)
      source_frames.push_back(i);
  else
    source_frames = parseRangeList<uint>(range);

  vector<uint> frames = trimFrames(source_frames, cutoff);
  if (frames.size() != source_frames.size())
    cout << "# WARNING- truncated last " << source_frames.size() - frames.size() << " frames\n";

  boost::tuple<vecGroup, vecUint> result = pickFiducials(subset, traj, frames, cutoff);
  cout << "# n\tref\n";
  vecGroup fiducials = boost::get<0>(result);
  vecUint id = boost::get<1>(result);

  for (uint i=0; i<id.size(); ++i)
    cout << i << "\t" << id[i] << "\n";


  DCDWriter dcd(outname + ".dcd", fiducials, hdr);

  PDB pdb = PDB::fromAtomicGroup(fiducials[0]);
  pdb.renumber();
  pdb.remarks().add(hdr);
  string fname = outname + ".pdb";
  ofstream ofs(fname.c_str());
  ofs << pdb;

}
