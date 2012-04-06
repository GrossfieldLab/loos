/*
  
  Perform a block-overlap in comparison to a full PCA

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

#include "ConvergenceOptions.hpp"
#include "bcomlib.hpp"


using namespace std;
using namespace loos;
using namespace Convergence;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;


typedef vector<AtomicGroup>                               vGroup;
typedef boost::tuple<RealMatrix, RealMatrix, RealMatrix>  SVDResult;





// Configuration

const bool length_normalize = true;
uint nsteps = 25;



// Global options
bool local_average;
bool use_zscore;
uint ntries;
vector<uint> blocksizes;
uint seed;
string gold_standard_trajectory_name;

string fullHelpMessage() {

  string s = 
    "\n"
    "SYNOPSIS\n"
    "\n"
    "Perform a block-overlap in comparison to a full PCA\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "This tool reports on how well a small \"block\" of a trajectory samples\n"
    "the subspace explored by the full simulation using principal component\n"
    "analysis.  It does this by computing the covariance overlap between a\n"
    "full simulation PCA and the PCA of increasingly longer contiguous \n"
    "blocks of that trajectory.\n "
    "\n"
    "See: Romo and Grossfield, J. Chem. Theor. Comput., 2011, 7, 2464-2472\n"
    "\n"
    "The output is a tab separated stream:\n"
    "n\tCoverlap\tVariance\tN_blocks\n"
    "\n"
    "\t   n     - current block size (nanoseconds)\n"
    "\tCoverlap - covariance overlap between block and full PCA\n"
    "\tVariance - variance in coverlap across all (N_blocks)\n"
    "\tN_blocks - number of blocks of a given length\n"
    "\n"
    "USAGE NOTES\n"
    "The --skip command is NOT used by this tool.\n"
    "\n"
    //
    "EXAMPLES\n"
    "bcom -s 'name==\"CA\"' --blocks 25:25:500 model.pdb traj.dcd > bcom_output\n"
    "\tCalculate the bcom of traj.dcd using a PCA of CA atoms. This\n"
    "\tis done for blocks in a range of 25 ns to 500 ns, with 25 ns\n"
    "\tintervals.  The result is written to the file bcom_output\n"
    "\n"
    "bcom -Z1 -s 'name==\"CA\"' --blocks 25:25:500 model.pdb traj.dcd > bcom_output\n"
    "\tSame as the example above, but outputs the block-averaged \n"
    "\tZ-score in the place of the block-averaged coverlap.\n"
    "\n"
    "bcom -s 'name==\"CA\"' --gold 'combined.dcd' model.pdb traj.dcd > bcom_output\n"
    "\tHere we make two changes.  First don't specify block sizes\n"
    "\tThis tells bcom to figure it out on its own.  In this case\n"
    "\tthe tool will run a max block size equal to half the trajectory.\n"
    "\tNext, we compare our block-averaged PCA results to a separate\n"
    "\ttrajectory called combined.dcd instead of the PCA of the full\n"
    "\ttraj.dcd.  As the name implies, combined.dcd may be a concatonation\n"
    "\tof several trajectories.  \n"
    "\t\tTo make such a concatoned trajectory see the tools\n"
    "\t\tmerge-traj and subsetter.\n"
    "\n"
    "* Visualization Notes *\n"
    "\tThe output should be plotted in the format X:Y:SQRT(Y-error)\n"
    "\twhere the colons separate the 1st 3 columns of the output.\n"
    "\tThis puts stdev error bars on the result\n"
    "\tIn GNUplot this would look like the following:\n"
    "\t   plot 'bcom_output' using 1:2:(sqrt(\\$3)) with errorlines\n"
    "\n";

  return(s);
}


// @cond TOOLS_INTERAL
class ToolOptions : public opts::OptionsPackage {
public:

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("blocks", po::value<string>(&blocks_spec), "Block sizes (MATLAB style range)")
      ("steps", po::value<uint>(&nsteps)->default_value(25), "Max number of blocks for auto-ranging")
      ("zscore,Z", po::value<bool>(&use_zscore)->default_value(false), "Use Z-score rather than covariance overlap")
      ("ntries,N", po::value<uint>(&ntries)->default_value(20), "Number of tries for Z-score")
      ("local", po::value<bool>(&local_average)->default_value(true), "Use local avg in block PCA rather than global")
      ("gold", po::value<string>(&gold_standard_trajectory_name)->default_value(""), "Use this trajectory for the gold-standard instead");

  }

  bool postConditions(po::variables_map& vm) {
    if (!blocks_spec.empty())
      blocksizes = parseRangeList<uint>(blocks_spec);

    return(true);
  }


  string print() const {
    ostringstream oss;
    oss << boost::format("blocks='%s', zscore=%d, ntries=%d, local=%d, gold='%s'")
      % blocks_spec
      % use_zscore
      % ntries
      % local_average
      % gold_standard_trajectory_name;
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






vGroup subgroup(const vGroup& A, const uint a, const uint b) {
  vGroup B;

  for (uint i=a; i<b; ++i)
    B.push_back(A[i]);

  return(B);
}


// Breaks the ensemble up into blocks and computes the PCA for each
// block and the statistics for the covariance overlaps...

template<class ExtractPolicy>
Datum blocker(const RealMatrix& Ua, const RealMatrix sa, vGroup& ensemble, const uint blocksize, ExtractPolicy& policy) {


  TimeSeries<double> coverlaps;

  for (uint i=0; i<ensemble.size() - blocksize; i += blocksize) {
    vGroup subset = subgroup(ensemble, i, i+blocksize);
    boost::tuple<RealMatrix, RealMatrix> pca_result = pca(subset, policy);
    RealMatrix s = boost::get<0>(pca_result);
    RealMatrix U = boost::get<1>(pca_result);

    // Scale the singular values by block-size
    if (length_normalize)
      for (uint j=0; j<s.rows(); ++j)
        s[j] /= blocksize;

    double val;
    if (use_zscore) {
      boost::tuple<double, double, double> result = zCovarianceOverlap(sa, Ua, s, U, ntries);
      val = boost::get<0>(result);
    } else
      val = covarianceOverlap(sa, Ua, s, U);

    coverlaps.push_back(val);
  }

  return( Datum(coverlaps.average(), coverlaps.variance(), coverlaps.size()) );
}




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicSelection* sopts = new opts::BasicSelection;
  opts::BasicTrajectory* tropts = new opts::BasicTrajectory;
  opts::BasicConvergence* copts = new opts::BasicConvergence;
  ToolOptions* topts = new ToolOptions;
  
  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(tropts).add(copts).add(topts);
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
 
  // First, align the input trajectory...
  boost::tuple<std::vector<XForm>, greal, int> ares = iterativeAlignment(ensemble);
  cout << "# Alignment converged to " << boost::get<1>(ares) << " in " << boost::get<2>(ares) << " iterations\n";
  cout << "# n\t" << (use_zscore ? "Z-score" : "Coverlap") << "\tVariance\tN_blocks\n";

  // Handle the gold-standard, either using the whole input traj, or the alternate traj...
  NoAlignPolicy policy;
  RealMatrix Us;
  RealMatrix UA;

  if (gold_standard_trajectory_name.empty()) {
    AtomicGroup avg = averageStructure(ensemble);
    policy = NoAlignPolicy(avg, local_average);
    boost::tuple<RealMatrix, RealMatrix> res = pca(ensemble, policy);

    Us = boost::get<0>(res);
    UA = boost::get<1>(res);

    if (length_normalize)
      for (uint i=0; i<Us.rows(); ++i)
        Us[i] /= traj->nframes();
  } else {
    // Must read in another trajectory, process it, and get the PCA
    pTraj gold = createTrajectory(gold_standard_trajectory_name, model);
    vector<AtomicGroup> gold_ensemble;
    readTrajectory(gold_ensemble, subset, gold);
    boost::tuple<vector<XForm>, greal, int> bres = iterativeAlignment(gold_ensemble);
    cout << "# Gold Alignment converged to " << boost::get<1>(bres) << " in " << boost::get<2>(bres) << " iterations\n";

    AtomicGroup avg = averageStructure(gold_ensemble);
    policy = NoAlignPolicy(avg, local_average);
    boost::tuple<RealMatrix, RealMatrix> res = pca(gold_ensemble, policy);

    Us = boost::get<0>(res);
    UA = boost::get<1>(res);


    if (length_normalize)
      for (uint i=0; i<Us.rows(); ++i)
        Us[i] /= gold->nframes();
  }



  // Now iterate over all requested block sizes

  // Provide user-feedback since this can be a slow computation
  PercentProgress watcher;
  ProgressCounter<PercentTrigger, EstimatingCounter> slayer(PercentTrigger(0.1), EstimatingCounter(blocksizes.size()));
  slayer.attach(&watcher);
  slayer.start();

  for (vector<uint>::iterator i = blocksizes.begin(); i != blocksizes.end(); ++i) {
    Datum result = blocker(UA, Us, ensemble, *i, policy);
    cout << *i << "\t" << result.avg_coverlap << "\t" << result.var_coverlap << "\t" << result.nblocks << endl;
    slayer.update();
  }

  slayer.finish();

}
