/*
  
  Perform block coverlap against self using z-score

*/



/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009-2011, Tod D. Romo
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
#include "bcomlib.hpp"

using namespace std;
using namespace loos;
using namespace Convergence;



typedef vector<AtomicGroup>                               vGroup;
typedef boost::tuple<RealMatrix, RealMatrix, RealMatrix>  SVDResult;


struct Datum {
  Datum(const double avg, const double var, const uint nblks) : average(avg),
                                                                variance(var),
                                                                nblocks(nblks) { }


  double average, variance;
  uint nblocks;
};



const bool length_normalize = true;




vGroup subgroup(const vGroup& A, const uint a, const uint b) {
  vGroup B;

  for (uint i=a; i<b; ++i)
    B.push_back(A[i]);

  return(B);
}



template<class ExtractPolicy>
Datum blocker(const uint n, vGroup& ensemble, const uint blocksize, ExtractPolicy& policy) {


  TimeSeries<double> zees;

  for (uint i=0; i<ensemble.size() - blocksize; i += blocksize) {
    vGroup subset = subgroup(ensemble, i, i+blocksize);
    boost::tuple<RealMatrix, RealMatrix> pca_result = pca(subset, policy);
    RealMatrix s = boost::get<0>(pca_result);
    RealMatrix U = boost::get<1>(pca_result);

    if (length_normalize)
      for (uint j=0; j<s.rows(); ++j)
        s[j] /= blocksize;

    boost::tuple<double, double, double> result = zCovarianceOverlap(s, U, s, U, n);
    double score = boost::get<0>(result);
    double cover = boost::get<1>(result);
    double dev = boost::get<2>(result);
    double avg = cover - dev * score;

    zees.push_back(avg / dev);
  }

  return( Datum(zees.average(), zees.variance(), zees.size()) );
}




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  int k=1;

  if (argc != 6) {
    cerr << "Usage- bcom model traj sel ntries blocks\n";
    exit(0);
  }

  AtomicGroup model = createSystem(argv[k++]);
  pTraj traj = createTrajectory(argv[k++], model);

  AtomicGroup subset = selectAtoms(model, argv[k++]);
  uint ntries = strtol(argv[k++], NULL, 10);
  vector<uint> blocksizes = parseRangeList<uint>(argv[k++]);

  vector<AtomicGroup> ensemble;
  readTrajectory(ensemble, subset, traj);
 
  // First, get the complete PCA result...
  boost::tuple<std::vector<XForm>, greal, int> ares = iterativeAlignment(ensemble);
  AtomicGroup avg = averageStructure(ensemble);
  NoAlignPolicy policy(avg, 1); // force local-avg

  cout << "# " << hdr << endl;
  cout << "# Config flags: length_normalize=" << length_normalize << endl;
  cout << "# Alignment converged to " << boost::get<1>(ares) << " in " << boost::get<2>(ares) << " iterations\n";
  cout << "# n\tZ-avg\tZ-var\tN_blocks\n";
  // Now iterate over all requested block sizes...
  
  PercentProgress watcher;
  ProgressCounter<PercentTrigger, EstimatingCounter> slayer(PercentTrigger(0.1), EstimatingCounter(blocksizes.size()));
  slayer.attach(&watcher);
  slayer.start();

  for (vector<uint>::iterator i = blocksizes.begin(); i != blocksizes.end(); ++i) {
    Datum result = blocker(ntries, ensemble, *i, policy);
    cout << boost::format("%d\t%10f\t%10f\t%d\n") % *i % result.average % result.variance % result.nblocks;
    slayer.update();
  }

  slayer.finish();

}
