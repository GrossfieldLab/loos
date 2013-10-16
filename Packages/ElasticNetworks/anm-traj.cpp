/*
  anm-traj

  (c) 2008,2013 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry

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

#include <limits>


#include "hessian.hpp"
#include "enm-lib.hpp"
#include "anm-lib.hpp"


using namespace std;
using namespace loos;
using namespace ENM;


namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;


// Globals...  Icky poo!

string prefix;

int verbosity;
bool debug;

string spring_desc;
string bound_spring_desc;

string fullHelpMessage() {

  string s = 
    "\n"                                                                                    
    "SYNOPSIS\n"                                                                      
    "\n"
    "ANM-based trajectory analysis (modeled after Hall, et al, JACS 129:11394 (2007))\n"
    "\n"                                                                                      
    "DESCRIPTION\n"                                                                      
    "\n"
    "Computes the anisotropic network model for each frame in a trajectory.\n"
    "The smallest non-zero eigenvalue is written to a matrix.  The corresponding\n"
    "eigenvector is also written as a column in another matrix.\n"
    "\n"
    "The following output files are created (using the optional prefix):\n"
    "\tgnm_traj_s.asc  - Smallest eigenvalue (magnitude of lowest frequency mode)\n"
    "\t                  First column is timestep, second column is the magnitude.\n"
    "\tgnm_traj_U.asc  - Corresponding eigenvectors.  Each column is an eigenvector\n"
    "\t                  and is paired with the corresponding row in the eigenvalue\n"
    "\t                  matrix file.\n"
    "\n"
    "\n"
    "* Spring Constant Control *\n"
    "Contacts between beads in an ANM are connected by a single potential\n"
    "which is described by a hookean spring.  The stiffness of each connection\n"
    "can be modified using various definitions of the spring constant.\n"
    "The spring constant used is controlled by the --spring option.\n"
    "If only the name for the spring function is given, then the default\n"
    "parameters are used.  Alternatively, the name may include a\n"
    "comma-separated list of parameters to be passed to the spring\n"
    "function, i.e. --spring=distance,15.0\n\n"
    "Available spring functions:\n";

  vector<string> names = springNames();
  for (vector<string>::const_iterator i = names.begin(); i != names.end(); ++i)
    s = s + "\t" + *i + "\n";
  
  s += 
    "\n\n"
    "* Adding \"Connectivity\" *\n"
    "ANM also supports construction of spring connections based on\n"
    "pseudo-connectivity.  This allows beads neighboring in sequence\n"
    "to be connected by a separate \"bound\" spring, chosen using the\n"
    "--bound option.  In this case the other or \"non-bound\" spring is\n"
    "chosen with the --spring option.\n"
    "\n"
    "\n\n"
    "EXAMPLES\n\n"
    "anm-traj --prefix b2ar b2ar.pdb b2ar.dcd\n"
    "\tCompute the ANM for all alpha-carbons in b2ar.  The output files are\n"
    "\tb2ar_s.asc (eigenvalues) and b2ar_U.asc (eigenvectors).\n"
    "\n"
    "anm-traj --selection 'resid >= 10 && resid <= 50 && name == \"CA\"' foo.pdb foo.dcd\n"
    "\tCompute the ANM for residues #10 through #50 with a 15 Angstrom cutoff\n"
    "\ti.e. construct contacts using only the CA's that are within 15 Angstroms\n"
    "\tThe model is in foo.pdb and the trajectory is stored in foo.dcd.  Output files\n"
    "\tcreated are anm_traj_s.asc (eigenvalues) and anm_traj_U.asc (eigenvectors).\n"
    "\n"
    "anm -S=exponential,-1.3 foo.pdb foo.dcd\n"
    "\tCompute an ANM using an spring function where the magnitude of\n"
    "\tthe connection decays exponentially with distance at a rate of\n"
    "\texp(-1.3*r) where r is the distance between contacts.  Note:\n"
    "\tin this case all beads are connected - which can eliminate\n"
    "\tan error in the numeric eigendecomposition.\n"
    "\n"
    "anm -b=constant,100 -S=exponential,-1.3 foo.pdb foo.dcd\n"
    "\tSimilar to the example above, but using connectivity.  Here\n"
    "\tresidues that are adjacent in sequence are connected by\n"
    "\tsprings with a constant stiffness of \"100\" and all other\n"
    "\tresidues are connected by springs that decay exponentially\n"
    "\twith distance\n"
    "\n"
    "NOTES\n"
    "- The default selection (if none is specified) is to pick CA's\n"
    "- The output is ASCII format suitable for use with Matlab/Octave/Gnuplot\n"
    "- Verbsity setting of 1 will give progress updates\n"
    "\n"
    "SEE ALSO\n"
    "\n"
    "gnm, gnm-traj, anm\n"
    "\n";

  return(s);
}




// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("debug", po::value<bool>(&debug)->default_value(false), "Turn on debugging (output intermediate matrices)")
      ("spring", po::value<string>(&spring_desc)->default_value("distance"),"Spring function to use")
      ("bound", po::value<string>(&bound_spring_desc), "Bound spring")
      ("vectors", po::value<bool>(&vectors)->default_value(false), "Write out matrix of first eigenvectors");
    
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("debug=%d, spring='%s', bound='%s', vectors=%d") % debug % spring_desc % bound_spring_desc % vectors;
    return(oss.str());
  }

  bool vectors;
};


class FastANM : public ElasticNetworkModel {
public:
  FastANM(SuperBlock* b) : ElasticNetworkModel(b) { prefix_ = "anm"; }

  void solve() {

    if (verbosity_ > 2)
      std::cerr << "Building hessian...\n";
    buildHessian();
    if (debugging_)
      loos::writeAsciiMatrix(prefix_ + "_H.asc", hessian_, meta_, false);

    loos::Timer<> t;
    if (verbosity_ > 1)
      std::cerr << "Computing Decomp of hessian...\n";
    t.start();

    eigenvals_ = eigenDecomp(hessian_);
    eigenvecs_ = hessian_;
    
    t.stop();
    if (verbosity_ > 1)
      std::cerr << "Decomp took " << loos::timeAsString(t.elapsed()) << std::endl;

  }


private:

};



// @endcond



loos::Math::Matrix<int> buildConnectivity(const AtomicGroup& model) {
  uint n = model.size();
  loos::Math::Matrix<int> M(n, n);
  
  for (uint j=0; j<n-1; ++j)
    for (uint i=j; i<n; ++i)
      if (i == j)
        M(j, i) = 1;
      else {
        M(j, i) = model[j]->isBoundTo(model[i]);
        M(i, j) = M(j, i);
      }
  
  return(M);
}


DoubleMatrix dotProduct(const DoubleMatrix& A) 
{
  uint m = A.rows();
  uint n = A.cols();
  
  DoubleMatrix D(n, n);

  for (uint j=0; j<n-1; ++j)
    for (uint i=j+1; i<n; ++i) {
      double s = 0.0;
      for (uint k=0; k<m; ++k)
  	s += A(k, j) * A(k, i);
      D(j, i) = D(i, j) = abs(s);
    }
  

  for (uint i=0; i<n; ++i)
    D(i, i) = 1.0;

  return(D);
}




int main(int argc, char *argv[]) {

  string header = invocationHeader(argc, argv);
  
  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::OutputPrefix* propts = new opts::OutputPrefix("anm_traj");
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

  verbosity = bopts->verbosity;
  prefix = propts->prefix;


  if (verbosity > 0)
    cerr << boost::format("Selected %d atoms from %s\n") % subset.size() % tropts->model_name;

  // Determine which kind of scaling to apply to the Hessian...
  vector<SpringFunction*> springs;
  SpringFunction* spring = 0;
  spring = springFactory(spring_desc);
  springs.push_back(spring);

  vector<SuperBlock*> blocks;
  SuperBlock* blocker = new SuperBlock(spring, subset);
  blocks.push_back(blocker);


  // Handle Decoration (if necessary)
  if (!bound_spring_desc.empty()) {
    if (! model.hasBonds()) {
      cerr << "Error- cannot use bound springs unless the model has connectivity\n";
      exit(-10);
    }
    loos::Math::Matrix<int> M = buildConnectivity(subset);
    SpringFunction* bound_spring = springFactory(bound_spring_desc);
    springs.push_back(bound_spring);

    BoundSuperBlock* decorator = new BoundSuperBlock(blocker, bound_spring, M);
    blocks.push_back(decorator);

    blocker = decorator;
  }


  FastANM anm(blocker);
  anm.debugging(debug);
  anm.prefix(prefix);
  anm.meta(header);
  anm.verbosity(verbosity);

  uint nframes = traj->nframes() - tropts->skip;
  uint natoms = subset.size();
  DoubleMatrix singvals(nframes, 3);
  DoubleMatrix singvecs(natoms, nframes);

  uint t = tropts->skip;
  uint k = 0;

  PercentProgressWithTime watcher;
  ProgressCounter<PercentTrigger, EstimatingCounter> slayer(PercentTrigger(0.1), EstimatingCounter(nframes - tropts->skip));
  slayer.attach(&watcher);
  if (verbosity)
    slayer.start();


  while (traj->readFrame()) {

    traj->updateGroupCoords(subset);
    anm.solve();

    DoubleMatrix s = anm.eigenvalues();
    singvals(k, 0) = t++;
    singvals(k, 1) = s[6];
    singvals(k, 2) = s[7];
    
    DoubleMatrix U = anm.eigenvectors();
    for (uint i=0; i<natoms; ++i)
      singvecs(i, k) = U(i, 6);

    if (verbosity)
      slayer.update();
    
    ++k;
  }

  writeAsciiMatrix(prefix + "_s.asc", singvals, header);
  if (topts->vectors)
    writeAsciiMatrix(prefix + "_U.asc", singvecs, header);

  DoubleMatrix D = dotProduct(singvecs);
  writeAsciiMatrix(prefix + "_D.asc", D, header);

  if (verbosity)
    slayer.finish();

  for (vector<SuperBlock*>::iterator i = blocks.begin(); i != blocks.end(); ++i)
    delete *i;
  for (vector<SpringFunction*>::iterator i = springs.begin(); i != springs.end(); ++i)
    delete *i;
  
}
