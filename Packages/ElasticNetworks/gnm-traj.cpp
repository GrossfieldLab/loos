/*
  gnm-traj

  Calculates a time-series of the first eigenvalue from a GNM calculated for each
  frame of a trajectory.

  See,
    Hall, B. A., Kaye, S. L., Pang, A., Perera, R. & Biggin, P. C. Characterization of protein conformational states by normal-mode frequencies. J Am Chem Soc 129, 11394–11401 (2007).

*/

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008,2013 Tod D. Romo
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


using namespace std;
using namespace loos;
namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;



string fullHelpMessage() {
  string msg = 
    "\n"
    "\n"
    "SYNOPSIS\n"
    "\n"
    "GNM-based trajectory analysis (see Hall, et al, JACS 129:11394 (2007))\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "Computes the gaussian network model for each frame in a trajectory.\n"
    "The smallest non-zero eigenvalue is written to a matrix.  The dot product\n"
    "of the corresponding eigenvector for each frame against every other frame\n"
    "is also written out.  The original eigenvectors may be optionally written as well.\n"
    "\n"
    "The following output files are created (using the optional prefix):\n"
    "\tgnm_traj_s.asc  - Smallest eigenvalue (magnitude of lowest frequency mode)\n"
    "\t                  First column is timestep, second column is the magnitude.\n"
    "\tgnm_traj_D.asc  - Matrix of dot products of corresponding eigenvectors.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "gnm-traj -v1 -pfoo -s 'resid >= 10 && resid <= 50 && name == \"CA\"' --cutoff 10.0 model.pdb traj.dcd\n"
    "\tPerform a GNM-analysis using model.pdb as the model and traj.dcd as the trajectory,\n"
    "\tfor residues #10 through #50 with a 10 Angstrom cutoff using only the C-alphas.\n"
    "\tWrites output files to foo_s.asc and foo_U.asc\n"
    "\tTiming and progress information will be written to the screen.\n"
    "\t\n"
    "NOTES\n"
    "- The default selection (if none is specified) is to pick CA's\n"
    "- The output is ASCII format suitable for use with Matlab/Octave/Gnuplot\n"
    "- Verbosity setting of 1 will give progress updates\n"
    "\n"
    "SEE ALSO\n"
    "\n"
    "gnm, anm, anm-traj\n"
    "\n";

  return(msg);
  
}


// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("cutoff", po::value<double>(&cutoff)->default_value(7.0), "Distance cutoff")
      ("vectors", po::value<bool>(&vectors)->default_value(false), "Write out matrix of first eigenvectors");
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("cutoff='%f',vectors=%d") % cutoff % vectors;
    return(oss.str());
  }

  double cutoff;
  bool vectors;
  
};


// @endcond


// This is the Kirchoff normalization constant (see Bahar, Atilgan,
// and Erman.  Folding & Design 2:173)
double normalization = 1.0;


DoubleMatrix kirchoff(AtomicGroup& group, const double cutoff) {
  int n = group.size();
  DoubleMatrix M(n, n);
  double r2 = cutoff * cutoff;


  for (int j=1; j<n; j++)
    for (int i=0; i<j; i++)
      if (group[i]->coords().distance2(group[j]->coords()) <= r2)
        M(i, j) = M(j, i) = -normalization;

  for (int j=0; j<n; j++) {
    double sum = 0;
    for (int i=0; i<n; i++) {
      if (i == j)
        continue;
      sum += M(j, i);
    }
    M(j, j) = -sum;
  }

  return(M);
}



DoubleMatrix dotProduct(const DoubleMatrix& A) 
{
  DoubleMatrix D = MMMultiply(A, A, true, false);
  for (ulong i=0; i<D.size(); ++i)
    D[i] = abs(D[i]);

  return(D);
}



int main(int argc, char *argv[]) 
{
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::OutputPrefix* propts = new opts::OutputPrefix("gnm_traj");
  opts::BasicSelection* sopts = new opts::BasicSelection("name == 'CA'");
  opts::BasicTrajectory* tropts = new opts::BasicTrajectory;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(propts).add(sopts).add(tropts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  AtomicGroup model = tropts->model;
  AtomicGroup subset = selectAtoms(model, sopts->selection);
  pTraj traj = tropts->trajectory;
  string prefix = propts->prefix;
  double cutoff = topts->cutoff;
  uint verbosity = bopts->verbosity;
  
  uint t = tropts->skip;
  uint k = 0;
  uint n = subset.size();
  DoubleMatrix svals(traj->nframes(), 3);
  DoubleMatrix vecs(n, traj->nframes());

  PercentProgressWithTime watcher;
  ProgressCounter<PercentTrigger, EstimatingCounter> slayer(PercentTrigger(0.1), EstimatingCounter(traj->nframes() - tropts->skip));
  slayer.attach(&watcher);
  if (verbosity)
    slayer.start();


  while (traj->readFrame()) {
    
    traj->updateGroupCoords(subset);
    DoubleMatrix K = kirchoff(subset, cutoff);
    DoubleMatrix S = eigenDecomp(K);
    
    svals(k, 0) = t++;
    svals(k, 1) = S[1];
    svals(k, 2) = S[2];

    for (uint i=0; i<n; ++i)
      vecs(i, k) = K(i, 1);
    
    ++k;

    if (verbosity)
      slayer.update();
  }

  writeAsciiMatrix(prefix + "_s.asc", svals, hdr);
  if (topts->vectors)
    writeAsciiMatrix(prefix + "_U.asc", vecs, hdr);

  DoubleMatrix D = dotProduct(vecs);
  writeAsciiMatrix(prefix + "_D.asc", D, hdr);

  if (verbosity)
    slayer.finish();

}


