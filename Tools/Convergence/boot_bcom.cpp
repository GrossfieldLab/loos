
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
#include <boost/program_options.hpp>
#include "bcomlib.hpp"

using namespace std;
using namespace loos;
using namespace Convergence;
namespace po = boost::program_options;


const bool debug = false;


typedef vector<AtomicGroup>                               vGroup;
typedef boost::tuple<RealMatrix, RealMatrix, RealMatrix>  SVDResult;

// Convenience structure for aggregating results
struct Datum {
  Datum(const double avg, const double var, const uint nblks) : avg_coverlap(avg),
                                                                var_coverlap(var),
                                                                nblocks(nblks) { }


  double avg_coverlap;
  double var_coverlap;
  uint nblocks;
};



// Configuration

const bool length_normalize = true;
const uint nsteps = 25;


// Global options
string model_name, traj_name, selection;
vector<uint> blocksizes;
bool local_average;
uint nreps;
uint seed;


void parseOptions(int argc, char *argv[]) {

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("seed,s", po::value<uint>(), "Seed for RNG")
      ("blocks,b", po::value<string>(), "Block sizes (MATLAB style range)")
      ("replicates,r", po::value<uint>(&nreps)->default_value(20), "Number of replicates for bootstrap")
      ("local,l", po::value<bool>(&local_average)->default_value(true), "Use local avg in block PCA rather than global");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "model")
      ("traj", po::value<string>(&traj_name), "trajectory")
      ("selection", po::value<string>(&selection), "selection");
    

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("traj", 1);
    p.add("selection", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("traj") && vm.count("selection"))) {
      cerr << "Usage- " << argv[0] << " [options] model trajectory selection >output\n";
      cerr << generic;
      exit(-1);
    }

    if (vm.count("blocks")) {
      string s = vm["blocks"].as<string>();
      blocksizes = parseRangeList<uint>(s);
    }

    if (vm.count("seed")) {
      seed = vm["seed"].as<uint>();
      rng_singleton().seed(seed);
    } else
      seed = randomSeedRNG();

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}



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

  parseOptions(argc, argv);

  cout << "# " << hdr << endl;
  cout << "# replicates = " << nreps << ", local_average = " << local_average << endl;
  cout << "# length_normalize = " << length_normalize << endl;

  AtomicGroup model = createSystem(model_name);
  pTraj traj = createTrajectory(traj_name, model);

  AtomicGroup subset = selectAtoms(model, selection);

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
