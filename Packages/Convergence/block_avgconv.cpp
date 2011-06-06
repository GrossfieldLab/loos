/*
  block_avgconv.cpp
  (c) 2010 Tod D. Romo, Grossfield Lab, URMC

  Convergence of average via block averaging

  Usage- block_avgconv model traj selection range
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
#include <boost/format.hpp>

using namespace loos;
using namespace std;


const uint default_starting_number_of_blocks = 500;
const double default_fraction_of_trajectory = 0.25;    


AtomicGroup averageSelectedSubset(const vector<AtomicGroup>& ensemble, const vector<uint>& indices) {
  AtomicGroup avg = ensemble[0].copy();
  for (AtomicGroup::iterator i = avg.begin(); i != avg.end(); ++i)
    (*i)->coords() = GCoord(0,0,0);

  uint n = avg.size();
  for (vector<uint>::const_iterator j = indices.begin(); j != indices.end(); ++j)
    for (uint i=0; i<n; ++i)
      avg[i]->coords() += ensemble[*j][i]->coords();

  for (AtomicGroup::iterator i = avg.begin(); i != avg.end(); ++i)
    (*i)->coords() /= indices.size();

  return(avg);
}




int main(int argc, char *argv[]) {
  
  if (argc < 4 || argc > 6) {
    cerr << "Usage- block_avgconv model traj sel [range [1 = do not align trajectory]]\n";
    exit(0);
  }

  string hdr = invocationHeader(argc, argv);

  int k = 1;
  AtomicGroup model = createSystem(argv[k++]);
  pTraj traj = createTrajectory(argv[k++], model);
  AtomicGroup subset = selectAtoms(model, argv[k++]);
  
  vector<uint> sizes;
  bool do_align = true;

  if (argc == 4) {
    uint step = traj->nframes() / default_starting_number_of_blocks;
    for (uint i=step; i<traj->nframes() * default_fraction_of_trajectory; i += step)
      sizes.push_back(i);
  } else {
    sizes = parseRangeList<uint>(argv[k++]);
    if (argc == 6)
      do_align = (argv[6][0] != '1');
  }


  cout << "# " << hdr << endl;
  cout << "# n\tavg\tvar\tblocks\tstderr\n";

  vector<AtomicGroup> ensemble;
  cerr << "Reading trajectory...\n";
  readTrajectory(ensemble, subset, traj);

  if (do_align) {
    cerr << "Aligning trajectory...\n";
    boost::tuple<vector<XForm>, greal, int> result = iterativeAlignment(ensemble);
  } else
    cerr << "Trajectory is already aligned!\n";

  cerr << "Processing- ";
  for (uint block = 0; block < sizes.size(); ++block) {
    if (block % 50)
      cerr << ".";

    uint blocksize = sizes[block];

    vector<AtomicGroup> averages;
    for (uint i=0; i<ensemble.size() - blocksize; i += blocksize) {

      vector<uint> indices(blocksize);
      for (uint j=0; j<blocksize; ++j)
        indices[j] = i+j;
      
      averages.push_back(averageSelectedSubset(ensemble, indices));
    }
    
    TimeSeries<double> rmsds;
    for (uint j=0; j<averages.size() - 1; ++j)
      for (uint i=j+1; i<averages.size(); ++i) {
        AtomicGroup left = averages[j];
        AtomicGroup right = averages[i];
        left.alignOnto(right);
        rmsds.push_back(left.rmsd(right));
      }
        

    double v = rmsds.variance();
    uint n = averages.size();
    cout << boost::format("%d\t%f\t%f\t%d\t%f\n") % blocksize % rmsds.average() % v % n % sqrt(v/n);
  }
  
  cerr << "\nDone!\n";
}
