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



int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  if (argc < 7 || argc > 8) {
    cerr << "Usage - fidpick model trajectory range|all selection output-name cutoff [seed]\n";
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

  vector<uint> frames;
  if (range == "all")
    for (uint i=0; i<traj->nframes(); ++i)
      frames.push_back(i);
  else
    frames = parseRangeList<uint>(range);

  boost::tuple<vecInt, vecUint, vecGroup, vecDouble> result = assignFrames(subset, traj, frames, cutoff);

  cout << "# " << hdr << endl;
  cout << "# seed = " << seed << endl;
  cout << boost::format("# i %-10s %-10s %-10s\n") % "Ref" % "Bin size" % "Radius";
  vecInt assignments = boost::get<0>(result);
  vecUint refs = boost::get<1>(result);
  vecGroup fiducials = boost::get<2>(result);
  vecDouble radii = boost::get<3>(result);

  vecUint hist = histogramBins(assignments);
  for (uint i=0; i<refs.size(); ++i)
    cout << boost::format("%-3d %-10d %-10d %10f\n") % i % refs[i] % hist[i] % radii[i];


  DCDWriter dcd(outname + ".dcd", fiducials, hdr);

  PDB pdb = PDB::fromAtomicGroup(fiducials[0]);
  pdb.renumber();
  pdb.remarks().add(hdr);
  string fname = outname + ".pdb";
  ofstream ofs(fname.c_str());
  ofs << pdb;

  fname = outname + "_assignments.asc";
  ofstream fass(fname.c_str());
  fass << "# " << hdr << endl;
  copy(assignments.begin(), assignments.end(), ostream_iterator<int>(fass, "\n"));

}
