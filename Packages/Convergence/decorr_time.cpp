/*
  Structural Histogram IID analysis a la Lyman & Zuckerman J Phys Chem
  B (2007) 111:12876-12882
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

#include "ConvergenceOptions.hpp"
#include "fid-lib.hpp"

using namespace std;
using namespace loos;


namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

// * GLOBALS *


// Turning on generates some internal debugging information
const bool debugging = false;
uint verbosity;

uint nreps = 5;
double frac;

vecUint trange;
vecUint nrange;
vecUint indices;

// @cond TOOLS_INTERNAL

string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\tCompute decorrelation times for a simulation\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tThis tool implements the decorrelation time method described in Lyman and Zuckerman,\n"
    "J Phys Chem B (2007) 111:12876-12882.  Bins for the structural histogram used are selected\n"
    "using a uniform probability, set with the --frac option (the default is 0.05 for 20 bins).\n"
    "The range of sample sizes (n, in figure 2) is given by the --nrange option, which takes\n"
    "either a comma-separated list of sizes or a matlab/octave-style range.  Finally, the\n"
    "required t-range is also a matlab/octave-style range (or comma-separated list) that\n"
    "gives the sample step-size (t in figure 2).  This is not to be confused with a range in\n"
    "frames of the trajectory.  However, the notion of \"time\" is dictated by the sampling\n"
    "rate of your trajectory, and is specified in terms of frames.  For example, if your\n"
    "trajectory has 1 frame/ns, then the t-range is specified in ns.  If your trajectory\n"
    "has one frame every 100 ps, then the t-range is specified in 100 ps units (i.e. frames).\n"
    "This whole procedure is repeated multiple times for each sample size, n.  The number of\n"
    "repeats is given by the --reps option (default of 5).\n"
    "\tThe output is an ASCII matrix where the first column is the step-size t. Each subsequent\n"
    "set of two-columns corresponds to a different sample size or n-value.  The first column\n"
    "in each set is the scaled variance (eq 3), averaged over each replica.  The second column\n"
    "is the standard error in the scaled variance.  This data can be plotted, e.g. figure 3.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\tdecorr_time --selection '!hydrogen' model.pdb simulation.dcd 5:5:100 >decorr.asc\n"
    "A decorrelation time calculation using the default sample sizes of 2, 4, and 10 and\n"
    "the default bin-size of 20 (probability 0.05).  The calculation is repeated the default\n"
    "of 5 times for each sample-size.  Only non-hydrogen atoms are used.  And the range in t\n"
    "is 5 to 100 at every 5 frames (assuming a 1 frame/ns trajectory, then 5 to 100 ns every 5 ns).\n"
    "\n"
    "\tdecorr_time --selection 'name == \"CA\"' --nrange 2,3,4 --frac 0.1 model.pdb simulation.dcd >decorr.asc 10:10:250\n"
    "Here, only alpha-carbons are used.  Sample sizes are 2, 3, and 4 and there are 10 bins.  The\n"
    "t-range here is 10 to 250 at every 10 frames.\n"
    "\n"
    "NOTES\n"
    "\tSimulations that are not well-converge may have difficulty with the default sample-size\n"
    "range.  Try a smaller range (i.e. 2,3,4) as well as large bins.\n"
    "SEE ALSO\n"
    "\tufidpick, effsize.pl, neff, assign_frames, hierarchy\n";

  return(msg);
}





class ToolOptions : public opts::OptionsPackage {
public:

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("nrange", po::value<string>(&nrange_spec)->default_value("2,4,10"), "Range of N to use")
      ("frac", po::value<double>(&frac)->default_value(0.05), "Bin fraction")
      ("reps", po::value<uint>(&nreps)->default_value(5), "# of repetitions to use for each N");
  }

  bool postConditions(po::variables_map& vm) {
    nrange = parseRangeList<uint>(nrange_spec);
    return(true);
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("nrange='%s', frac=%f, reps=%f")
      % nrange_spec
      % frac
      % nreps;
    return(oss.str());
  }

  string nrange_spec;
};


// @endcond





// Despite the name, this histograms the assignments (i.e. frequency
// of the various fiducials in the trajectory)

vecDouble rebinFrames(const vecUint& assignments, const uint S, const vecUint& ensemble) {
  vecUint hist(S, 0);

  for (vecUint::const_iterator i = ensemble.begin(); i != ensemble.end(); ++i) {
    if (assignments[*i] >= S)
      throw(runtime_error("Bin index exceeds number of fiducials in rebinFrames()"));

    ++hist[assignments[*i]];
  }

  vecDouble pops(S);
  for (uint i=0; i<S; ++i)
    pops[i] = static_cast<double>(hist[i])/ensemble.size();

  return(pops);
}




vecDouble binVariances(const vecUint& assignments, const uint S, const uint n, const uint t) {
  vector<vecDouble> fik;
  
  uint curframe = 0;
  while (curframe < assignments.size()) {
    vecUint ensemble;
    for (uint i=0; i<n && curframe < assignments.size(); ++i, curframe += t)
      ensemble.push_back(curframe);
         
    // Incomplete sample, so ignore...
    if (ensemble.size() < n)
      break;

    vecDouble hist = rebinFrames(assignments, S, ensemble);
    fik.push_back(hist);
  }

  vecDouble means(S, 0.0);
  for (vector<vecDouble>::iterator j = fik.begin(); j != fik.end(); ++j)
    for (uint i=0; i<S; ++i)
      means[i] += (*j)[i];
  for (uint i=0; i<S; ++i)
    means[i] /= fik.size();

  if (debugging) {
    cerr << "Probe> means=";
    copy(means.begin(), means.end(), ostream_iterator<double>(cerr, ","));
    cerr << endl;
  }

  vecDouble vars(S, 0.0);

  for (vector<vecDouble>::iterator j = fik.begin(); j != fik.end(); ++j)
    for (uint i=0; i<S; ++i) {
      double d = (*j)[i] - means[i];
      vars[i] += d*d;
    }

  for (vecDouble::iterator i = vars.begin(); i != vars.end(); ++i)
    *i /= fik.size();

  if (debugging) {
    cerr << boost::format("Probe> chunks = %d\n") % fik.size();
    cerr << "Probe> vars=";
    copy(vars.begin(), vars.end(), ostream_iterator<double>(cerr, ","));
    cerr << endl;
  }

  return(vars);
}


double avgVariance(const vecDouble& vars) {
  double mean = 0.0;
  
  for (vecDouble::const_iterator i = vars.begin(); i != vars.end(); ++i)
    mean += *i;
  mean /= vars.size();

  return(mean);
}


double sigma(const vecUint& assignments, const uint S, const uint n, const uint t) {
  
  vecDouble vars = binVariances(assignments, S, n, t);
  double mean_vars = avgVariance(vars);
  
  double f = 1.0 / S;
  double N = static_cast<double>(assignments.size()) / t;
  double expected = (f*(1.0-f) / n) * (N-n)/(N-1.0);

  if (debugging)
    cerr << boost::format("Probe> f=%f, N=%f, expected=%f\n") % f % N % expected;
  return(mean_vars/expected);
}


DoubleMatrix statistics(const vector<DoubleMatrix>& V) {
  uint m = V[0].rows();
  uint n = V[0].cols();

  DoubleMatrix M(m, n);
  for (uint k=0; k<V.size(); ++k)
    for (ulong i=0; i<m*n; ++i)
      M[i] += V[k][i];
  for (long i=0; i<m*n; ++i)
    M[i] /= V.size();

  DoubleMatrix S(m, n);
  if (V.size() > 2) {
    for (uint k=0; k<V.size(); ++k)
      for (ulong i=0; i<m*n; ++i) {
        double d = V[k][i] - M[i];
        S[i] = d*d;
      }
    for (ulong i=0; i<m*n; ++i)
      S[i] = sqrt(S[i] / (V.size() - 1));
  }

  uint nn = (n-1) * 2 + 1;
  DoubleMatrix R(m, nn);
  for (uint j=0; j<m; ++j)
    for (uint i=0; i<n; ++i)
      if (i == 0)
        R(j, i) = M(j, i);
      else {
        R(j, 2*(i-1)+1) = M(j, i);
        R(j, 2*(i-1)+2) = S(j, i);
      }
  
  return(R);
}



int main(int argc, char *argv[]) {

  string hdr = invocationHeader(argc, argv);
  opts::BasicOptions *bopts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicSelection *sopts = new opts::BasicSelection;
  opts::TrajectoryWithFrameIndices *tropts = new opts::TrajectoryWithFrameIndices;
  opts::BasicConvergence *copts = new opts::BasicConvergence;
  ToolOptions* topts = new ToolOptions;
  opts::RequiredArguments* ropts = new opts::RequiredArguments("trange", "T-range");
  
  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(tropts).add(copts).add(topts).add(ropts);

  if (!options.parse(argc, argv))
    exit(-1);

  verbosity = bopts->verbosity;

  cout << "# " << hdr << endl;
  cout << "# " << vectorAsStringWithCommas<string>(options.print()) << endl;


  AtomicGroup model = tropts->model;
  pTraj traj = tropts->trajectory;
  AtomicGroup subset = selectAtoms(model, sopts->selection);

  string trange_spec = ropts->value("trange");
  trange = parseRangeList<uint>(trange_spec);

  indices = assignTrajectoryFrames(traj, tropts->frame_index_spec, tropts->skip);

  vector<DoubleMatrix> results;
  for (uint k = 0; k<nreps; ++k) {
    if (verbosity > 0)
      cerr << "Replica #" << k << endl;

    boost::tuple<vecGroup, vecUint> fids = pickFiducials(subset, traj, indices, frac);
    vecGroup fiducials = boost::get<0>(fids);
    vecUint assignments = assignStructures(subset, traj, indices, fiducials);
    uint S = fiducials.size();
    
    DoubleMatrix M(trange.size(), nrange.size() + 1);
    for (uint i=0; i<nrange.size(); ++i)
      for (uint j=0; j<trange.size(); ++j)
        M(j, i+1) = sigma(assignments, S, nrange[i], trange[j]);
    
    for (uint j=0; j<trange.size(); ++j)
      M(j, 0) = trange[j];

    results.push_back(M);
  }

  DoubleMatrix M = statistics(results);
  writeAsciiMatrix(cout, M, "");
}

