/*
  avgconv

  Usage- avgconv model traj selection range
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
#include <boost/format.hpp>

using namespace std;
using namespace loos;


bool locally_optimal = false;


AtomicGroup calcAverage(const vector<AtomicGroup>& ensemble, const uint size) {

  vector<AtomicGroup> subsample(size);
  for (uint i=0; i<size; ++i)
    subsample[i] = ensemble[i];

  if (locally_optimal)
    (void)iterativeAlignment(subsample);

  AtomicGroup avg = averageStructure(subsample);
  return(avg);
}



int main(int argc, char *argv[]) {
  if (argc < 4 || argc > 6) {
    cout << "Usage- avgconv model traj selection [range [1 = local optimal avg]]\n";
    exit(-1);
  }

  string hdr = invocationHeader(argc, argv);
  cout << "# " << hdr << endl;
  cout << "# n\trmsd\n";

  int k = 1;
  AtomicGroup model = createSystem(argv[k++]);
  pTraj traj = createTrajectory(argv[k++], model);
  string sel = string(argv[k++]);
  
  vector<uint> blocks;
  if (argc == 4) {

    uint step = traj->nframes() / 100;
    for (uint i=step; i<traj->nframes(); i += step)
      blocks.push_back(i);

  } else {

    blocks = parseRangeList<uint>(argv[k++]);
    if (argc == 6)
      locally_optimal = (argv[6][0] == '1');
  }
  
  AtomicGroup subset = selectAtoms(model, sel);
  cout << boost::format("# Subset has %d atoms\n") % subset.size();
  vector<AtomicGroup> ensemble;
  readTrajectory(ensemble, subset, traj);
  cout << boost::format("# Trajectory has %d frames\n") % ensemble.size();

  cout << boost::format("# Blocks = %d\n") % blocks.size();

  if (!locally_optimal) {
    boost::tuple<vector<XForm>, greal, int> result = iterativeAlignment(ensemble);
    cout << boost::format("# Iterative alignment converged to RMSD of %g with %d iterations\n") % boost::get<1>(result) % boost::get<2>(result);
  }

  AtomicGroup preceding = calcAverage(ensemble, blocks[0]);
  for (vector<uint>::const_iterator ci = blocks.begin()+1; ci != blocks.end(); ++ci) {
    AtomicGroup avg = calcAverage(ensemble, *ci);
    avg.alignOnto(preceding);
    double rmsd = preceding.rmsd(avg);

    cout << boost::format("%d\t%f\n") % *ci % rmsd;

    preceding = avg;
  }

}

