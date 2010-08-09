/*
HI NICK!!!!
  anm

  (c) 2008 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry

  Computes the anisotropic network model for a structure.  It does
  this by building a hessian for the structure, then computing the SVD
  of it and the corresponding pseudo-inverse (ignoring the 6 lowest
  modes).

  Usage:
    anm [selection string] radius model-name output-prefix

  Examples:
    anm 'resid >= 10 && resid <= 50 && name == "CA"' 15.0 foo.pdb foo

    This creates the following files:
      foo_H.asc      == The hessian
      foo_U.asc      == Left singular vectors
      foo_s.asc      == Singular values
      foo_V.asc      == Right singular vectors
      foo_Hi.asc     == Pseudo-inverse of H

  Notes:
    o The default selection (if none is specified) is to pick CA's
    o The output is in ASCII format suitable for use with Matlab/Octave/Gnuplot
      
  
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008,2009 Tod D. Romo
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
//#include "oldhessian.hpp"

#include "hessian.hpp"
#include "enm-lib.hpp"
#include "anm-lib.hpp"

#include <limits>
#include <boost/format.hpp>
#include <boost/program_options.hpp>


using namespace std;
using namespace loos;
namespace po = boost::program_options;


typedef Math::Matrix<double, Math::ColMajor> Matrix;



// Globals...  Icky poo!

string selection;
string bsf, nbsf;

string model_name;
string prefix;
double cutoff;

// Turns on parameter-free mode a la Yang et al, PNAS (2009) 106:12347
bool parameter_free;

double power;
bool hca_method;

int verbosity;
bool debug;


double hca_constants[5];
bool user_defined_hca_constants(false);


void fullHelp() {

  cout <<
    "\n"
    "Computes the anisotropic network model for a structure.  It does\n"
    "this by building a hessian for the structure, then computing the SVD\n"
    "of it and the corresponding pseudo-inverse (ignoring the 6 lowest\n"
    "modes).\n"
    "\n"
    "This creates the following files:\n"
    "\tfoo_H.asc      == The hessian\n"
    "\tfoo_U.asc      == Left singular vectors\n"
    "\tfoo_s.asc      == Singular values\n"
    "\tfoo_V.asc      == Right singular vectors\n"
    "\tfoo_Hi.asc     == Pseudo-inverse of H\n"
    "\n"
    "\n"
    "* Spring Constant Control *\n\n"
    "Different methods for assigning the spring constants in the\n"
    "Hessian can be used.  For example, \"--free 1\" selects the\n"
    "\"parameter free\" method which can be combined with the \"--power\"\n"
    "option, which controls the exponent used (the default is -2).\n"
    "Note that setting the parameter-free method does not alter the\n"
    "cutoff radius used in building the Hessian, so you may want to\n"
    "set that to something very large (i.e. \"--cutoff 100\").\n"
    "Alternatively, the \"HCA\" method can be used via the \"--hca 1\"\n"
    "option.  The constants used in HCA can be set with the\n"
    "\"--hparams r_c,k1,k2,k3,k4\" option where the spring constant, k,\n"
    "is defined as,\n"
    "\tk = k1 * s - k2        if (s <= r_c)\n"
    "\tk = k3 * pow(s, -k4)   if (s > r_c)\n"
    "and s is the distance between the nodes.\n"
    "\n\n"
    "Examples:\n\n"
    "Compute the ANM for residues #10 through #50 with a 15 Angstrom cutoff\n"
    "\tanm 'resid >= 10 && resid <= 50 && name == \"CA\"' 15.0 foo.pdb foo\n";
}



void parseOptions(int argc, char *argv[]) {

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("verbosity,v", po::value<int>(&verbosity)->default_value(0), "Verbosity level")
      ("debug,d", po::value<bool>(&debug)->default_value(false), "Turn on debugging (output intermediate matrices)")
      ("selection,s", po::value<string>(&selection)->default_value("name == 'CA'"), "Which atoms to use for the network")
      ("free,f", po::value<bool>(&parameter_free)->default_value(false), "Use the parameter-free method rather than the cutoff")
      ("hca,h", po::value<bool>(&hca_method)->default_value(false), "Use the HCA distance scaling method")
      ("hparams,H", po::value<string>(), "Constants to use in HCA scaling (rcut, k1, k2, k3, k4)")
      ("power,P", po::value<double>(&power)->default_value(-2.0), "Scale to use for parameter-free")
      ("cutoff,c", po::value<double>(&cutoff)->default_value(15.0), "Cutoff distance for node contact")
      ("fullhelp", "More detailed help")
      ("bonded_function,b", po::value<string>(&bsf)->default_value("exponential"), "Which spring funtion should be used for bonded nodes")
      ("nonbonded_function,n", po::value<string>(&nbsf)->default_value("exponential"), "Which spring funtion should be used for NONbonded nodes");


    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("prefix", po::value<string>(&prefix), "Output prefix");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("prefix", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || vm.count("fullhelp") || !(vm.count("model") && vm.count("prefix"))) {
      cerr << "Usage- anm [options] model-name output-prefix\n";
      cerr << generic;
      if (vm.count("fullhelp"))
        fullHelp();
      exit(-1);
    }

    // Force the hessian to include all nodes...
    if (parameter_free)
      cutoff = numeric_limits<double>::max();

    // if (vm.count("bonded_function"))
    //   string bsf = vm["bonded_function"].as<string>();
    
    // if (vm.count("nonbonded_function"))
    //   string nbsf = vm["nonbonded_function"].as<string>();

    if (vm.count("hparams")) {
      string s = vm["hparams"].as<string>();
      int i = sscanf(s.c_str(), "%lf,%lf,%lf,%lf,%lf", hca_constants, hca_constants+1, hca_constants+2,
                     hca_constants+3, hca_constants+4);
      if (i != 5) {
        cerr << boost::format("Error - invalid conversion of HCA constants '%s'\n") % s;
        exit(-1);
      }
      user_defined_hca_constants = true;
    }


  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}


int main(int argc, char *argv[]) {

  
  string header = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  AtomicGroup model = createSystem(model_name);
  AtomicGroup subset = selectAtoms(model, selection);

  if (verbosity > 0)
    cerr << boost::format("Selected %d atoms from %s\n") % subset.size() % model_name;

  // Determine which kind of scaling to apply to the Hessian...
  SpringFunction* bound_spring = 0;
  SpringFunction* nonbound_spring = 0;

  if (!bsf.empty()) 
    bound_spring = springFactory(bsf);
  if (!nbsf.empty())
    nonbound_spring  = springFactory(nbsf);

  cerr << "\n made springs...\t" << subset.size();
//   Adding the connectivity map
  loos::Math::Matrix<int> connectivity_map(subset.size(), subset.size());
  if (subset.hasBonds()){
    cerr << "\n has bonds!";
    for (int j = 0; subset.size(); ++j){
      cerr << "\n current j:  " << j;
      if (subset[j].hasBonds()){
	vector<int> jbonds = subset[j]->getBonds();
	// jbonds now holds the indices for all atoms bonded to j
	for (int i = 0; i < jbonds.size(); ++i){
	  cerr << "\n current i:  " << i;
	  for (int k = 0; k < subset.size(); ++k) {
	    // cerr << "\t current k:  " << k;
	    if (jbonds[i] == subset[k]->id())
	      connectivity_map(j,k) = 1;
	    else
	      connectivity_map(j,k) = 0;
	  }
	}
      }
    }
  }
  cerr << "defined connectivity map \n";
  //   Impleminting the decorator
  SuperBlock* forbondedTerms = new SuperBlock(bound_spring, subset);
  BoundSuperBlock* forAllTerms = new BoundSuperBlock(forbondedTerms, nonbound_spring, connectivity_map);
  //loop over connectivity_map
  //if (j,i) = 1 than use bondedTerm if = 0 use nonbondedTerm
  //does this get rid of the need for DoubleMatrix H??


  /////i think that this should be a boundsuperblock!!!!
  ANM anm(forAllTerms);
  anm.prefix(prefix);
  anm.meta(header);
  anm.solve();
  
  writeAsciiMatrix(prefix + "_map.asc", connectivity_map, header, false);

  // Write out the LSVs (or eigenvectors)
  writeAsciiMatrix(prefix + "_U.asc", anm.eigenvectors(), header, false);
  writeAsciiMatrix(prefix + "_s.asc", anm.eigenvalues(), header, false);

  writeAsciiMatrix(prefix + "_Hi.asc", anm.inverseHessian(), header, false);

  delete forAllTerms;
  delete forbondedTerms;
}
