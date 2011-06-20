/*
  
  Cosine content for varying windows of a trajectory,
  based on:
    Hess, B.  "Convergence of sampling in protein simulations."
      Phys Rev E (2002) 65(3):031910

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

#include "bcomlib.hpp"

using namespace std;
using namespace loos;
using namespace Convergence;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;


typedef vector<AtomicGroup>                               vGroup;


// Convenience structure for aggregating results
struct Datum {
  Datum(const double avg, const double var, const uint nblks) : avg_cosine(avg),
                                                                var_cosine(var),
                                                                nblocks(nblks) { }


  double avg_cosine;
  double var_cosine;
  uint nblocks;
};



// Configuration

const uint nsteps = 50;



// Global options
bool local_average;
vector<uint> blocksizes;
string model_name, traj_name, selection;
uint principal_component;


// @cond TOOLS_INTERAL
class ToolOptions : public opts::OptionsPackage {
public:

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("pc", po::value<uint>(&principal_component)->default_value(0), "Which principal component to use")
      ("blocks", po::value<string>(&blocks_spec), "Block sizes (MATLAB style range)")
      ("local", po::value<bool>(&local_average)->default_value(true), "Use local avg in block PCA rather than global");

  }

  bool postConditions(po::variables_map& vm) {
    if (!blocks_spec.empty())
      blocksizes = parseRangeList<uint>(blocks_spec);

    return(true);
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("blocks='%s', local=%d, pc=%d")
      % blocks_spec
      % local_average
      % principal_component;
    return(oss.str());
  }

  string blocks_spec;
};
// @endcond







vGroup subgroup(const vGroup& A, const uint a, const uint b) {
  vGroup B;

  for (uint i=a; i<b; ++i)
    B.push_back(A[i]);

  return(B);
}



// Breaks the ensemble up into blocks and computes the RSV for each
// block and the statistics for the cosine content...

template<class ExtractPolicy>
Datum blocker(const uint pc, vGroup& ensemble, const uint blocksize, ExtractPolicy& policy) {


  TimeSeries<double> cosines;

  for (uint i=0; i<ensemble.size() - blocksize; i += blocksize) {
    vGroup subset = subgroup(ensemble, i, i+blocksize);
    RealMatrix V = rsv(subset, policy);

    double val = cosineContent(V, pc);
    cosines.push_back(val);
  }

  return( Datum(cosines.average(), cosines.variance(), cosines.size()) );
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
    uint step = n / nsteps;
    if (step < 1)
      step = 1;
    cout << "# Auto block-sizes - " << step << ":" << step << ":" << n-1 << endl;

    for (uint i = step; i < n; i += step)
      blocksizes.push_back(i);
  }


  AtomicGroup subset = selectAtoms(model, sopts->selection);


  vector<AtomicGroup> ensemble;
  readTrajectory(ensemble, subset, traj);
 
  // First, read in and align trajectory
  boost::tuple<std::vector<XForm>, greal, int> ares = iterativeAlignment(ensemble);
  AtomicGroup avg = averageStructure(ensemble);
  NoAlignPolicy policy(avg, local_average);


  // Now iterate over all requested block sizes

  // Provide user-feedback since this can be a slow computation
  PercentProgress watcher;
  ProgressCounter<PercentTrigger, EstimatingCounter> slayer(PercentTrigger(0.1), EstimatingCounter(blocksizes.size()));
  slayer.attach(&watcher);
  slayer.start();

  for (vector<uint>::iterator i = blocksizes.begin(); i != blocksizes.end(); ++i) {
    Datum result = blocker(principal_component, ensemble, *i, policy);
    cout << *i << "\t" << result.avg_cosine << "\t" << result.var_cosine << "\t" << result.nblocks << endl;
    slayer.update();
  }

  slayer.finish();

}
