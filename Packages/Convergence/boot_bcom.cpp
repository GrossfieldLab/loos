
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

namespace opts = loos::OptionsFramework;
namespace po = boost::program_options;


const bool debug = false;


typedef vector<AtomicGroup>                               vGroup;
typedef boost::tuple<RealMatrix, RealMatrix, RealMatrix>  SVDResult;


// Configuration

const bool length_normalize = true;
const uint nsteps = 25;


// Global options
vector<uint> blocksizes;
bool local_average;
uint nreps;
uint seed;

// @cond TOOLS_INTERAL
class ToolOptions : public opts::OptionsPackage {
public:

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("blocks", po::value<string>(&blocks_spec), "Block sizes (MATLAB style range)")
      ("reps", po::value<uint>(&nreps)->default_value(20), "Number of replicates for bootstrap")
      ("local", po::value<bool>(&local_average)->default_value(true), "Use local avg in block PCA rather than global")
      ("seed", po::value<uint>(&seed)->default_value(0), "Random number seed (0 = auto)");

  }

  bool postConditions(po::variables_map& vm) {
    if (!blocks_spec.empty())
      blocksizes = parseRangeList<uint>(blocks_spec);

    if (seed == 0)
      seed = randomSeedRNG();
    else
      rng_singleton().seed(seed);

    return(true);
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("blocks='%s', local=%d, reps=%d, seed=%d")
      % blocks_spec
      % local_average
      % nreps
      % seed;
    return(oss.str());
  }

  string blocks_spec;
};

// Convenience structure for aggregating results
struct Datum {
  Datum(const double avg, const double var, const uint nblks) : avg_coverlap(avg),
                                                                var_coverlap(var),
                                                                nblocks(nblks) { }


  double avg_coverlap;
  double var_coverlap;
  uint nblocks;
};
// @endcond


// Randomly pick frames
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


// Extract a subgroup of the vector<AtomicGroup> given the indices in picks...
vGroup subgroup(const vGroup& A, const vector<uint>& picks) {
  vGroup B;
  
  for (vector<uint>::const_iterator ci = picks.begin(); ci != picks.end(); ++ci)
    B.push_back(A[*ci]);

  return(B);
}



// Breaks the ensemble up into blocks and computes the PCA for each
// block and the statistics for the covariance overlaps...

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

  opts::BasicOptions* bopts = new opts::BasicOptions;
  opts::BasicSelection* sopts = new opts::BasicSelection;
  opts::BasicTrajectory* tropts = new opts::BasicTrajectory;
  ToolOptions* topts = new ToolOptions;
  
  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(tropts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  cout << "# " << hdr << endl;
  cout << "# " << vectorAsStringWithCommas<string>(options.print()) << endl;

  AtomicGroup model = tropts->model;
  pTraj traj = tropts->trajectory;
  if (tropts->skip)
    cerr << "Warning: --skip option ignored\n";


  if (blocksizes.empty()) {
    uint n = traj->nframes();
    uint half = n / 2;
    uint step = half / nsteps;
    if (step < 1)
      step = 1;
    cout << "# Auto block-sizes - " << step << ":" << step << ":" << half << endl;

    for (uint i = step; i <= half; i += step)
      blocksizes.push_back(i);
  }


  AtomicGroup subset = selectAtoms(model, sopts->selection);


  vector<AtomicGroup> ensemble;
  readTrajectory(ensemble, subset, traj);

  // First, get the complete PCA result...
  boost::tuple<std::vector<XForm>, greal, int> ares = iterativeAlignment(ensemble);
  AtomicGroup avg = averageStructure(ensemble);
  NoAlignPolicy policy(avg, local_average);
  boost::tuple<RealMatrix, RealMatrix> res = pca(ensemble, policy);

  RealMatrix Us = boost::get<0>(res);
  RealMatrix UA = boost::get<1>(res);


  if (length_normalize)
    for (uint i=0; i<Us.rows(); ++i)
      Us[i] /= traj->nframes();


  
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
