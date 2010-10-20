
/*
  Perform a bootstrap analysis of a trajectory
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
#include "bcomlib.hpp"

using namespace std;
using namespace loos;
using namespace Convergence;


const bool debug = false;


typedef vector<AtomicGroup>                               vGroup;
typedef boost::tuple<RealMatrix, RealMatrix, RealMatrix>  SVDResult;


struct Datum {
  Datum(const double avg, const double var, const uint nblks) : avg_coverlap(avg),
                                                                var_coverlap(var),
                                                                nblocks(nblks) { }


  double avg_coverlap;
  double var_coverlap;
  uint nblocks;
};




const bool length_normalize = true;



vector<uint> pickFrames(const uint nframes, const uint blocksize) {
  
  boost::uniform_int<uint> imap(0,nframes-1);
  boost::variate_generator< base_generator_type&, boost::uniform_int<uint> > rng(rng_singleton(), imap);
  vector<uint> picks;

  for (uint i=0; i<blocksize; ++i)
    picks.push_back(rng());

  return(picks);
}


void dumpPicks(const vector<uint>& picks) {
  cerr << "Picks:\n";
  for (vector<uint>::const_iterator ci = picks.begin(); ci != picks.end(); ++ci)
    cerr << "\t" << *ci << endl;
}


vGroup subgroup(const vGroup& A, const vector<uint>& picks) {
  vGroup B;
  
  for (vector<uint>::const_iterator ci = picks.begin(); ci != picks.end(); ++ci)
    B.push_back(A[*ci]);

  return(B);
}



template<class ExtractPolicy>
Datum blocker(const RealMatrix& Ua, const RealMatrix sa, const vGroup& ensemble, const uint blocksize, uint repeats, ExtractPolicy& policy) {


  
  TimeSeries<double> coverlaps;

  for (uint i=0; i<repeats; ++i) {
    vector<uint> picks = pickFrames(ensemble.size(), blocksize);
    
    if (debug) {
      cerr << "***Block " << blocksize << ", replica " << i << ", picks " << picks.size() << endl;
      dumpPicks(picks);
    }
    
    vGroup subset = subgroup(ensemble, picks);
    boost::tuple<RealMatrix, RealMatrix> pca_result = pca(subset, policy);
    RealMatrix s = boost::get<0>(pca_result);
    RealMatrix U = boost::get<1>(pca_result);

    if (length_normalize)
      for (uint j=0; j<s.rows(); ++j)
        s[j] /= blocksize;

    coverlaps.push_back(covarianceOverlap(sa, Ua, s, U));
  }

  return( Datum(coverlaps.average(), coverlaps.variance(), coverlaps.size()) );

}




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  int k=1;

  if (argc != 8) {
    cerr << "Usage- " << argv[0] << " model traj sel replicates [0|seed] [1=local avg|0=global avg] blocks\n";
    exit(0);
  }

  AtomicGroup model = createSystem(argv[k++]);
  pTraj traj = createTrajectory(argv[k++], model);

  AtomicGroup subset = selectAtoms(model, argv[k++]);
  int nreps = strtol(argv[k++],0, 10);
  uint seed = strtol(argv[k++],0, 10);
  if (seed == 0)
    seed = randomSeedRNG();
  else
    rng_singleton().seed(static_cast<ulong>(seed));

  int local_flag = atoi(argv[k++]);

  vector<uint> blocksizes = parseRangeList<uint>(argv[k++]);

  vector<AtomicGroup> ensemble;
  readTrajectory(ensemble, subset, traj);

  // First, get the complete PCA result...
  boost::tuple<std::vector<XForm>, greal, int> ares = iterativeAlignment(ensemble);
  AtomicGroup avg = averageStructure(ensemble);
  NoAlignPolicy policy(avg, local_flag);
  boost::tuple<RealMatrix, RealMatrix> res = pca(ensemble, policy);

  RealMatrix Us = boost::get<0>(res);
  RealMatrix UA = boost::get<1>(res);


  if (length_normalize)
    for (uint i=0; i<Us.rows(); ++i)
      Us[i] /= traj->nframes();


  
  cout << "# " << hdr << endl;
  cout << "# Config flags: length_normalize=" << length_normalize << endl;
  cout << "# Alignment converged to " << boost::get<1>(ares) << " in " << boost::get<2>(ares) << " iterations\n";
  cout << "# seed = " << seed << endl;
  cout << "# n\tCoverlap\tVariance\tN_blocks\n";
  // Now iterate over all requested block sizes...

  PercentProgress watcher;
  ProgressCounter<PercentTrigger, EstimatingCounter> slayer(PercentTrigger(0.1), EstimatingCounter(blocksizes.size()));
  slayer.attach(&watcher);
  slayer.start();


  for (vector<uint>::iterator i = blocksizes.begin(); i != blocksizes.end(); ++i) {
    Datum result = blocker(UA, Us, ensemble, *i, nreps, policy);
    cout << *i << "\t" << result.avg_coverlap << "\t" << result.var_coverlap << "\t" << result.nblocks << endl;
    slayer.update();
  }

  slayer.finish();
}
