/*
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

#include <limits>
#include <boost/format.hpp>
#include <boost/program_options.hpp>


#include "hessian.hpp"
#include "enm-lib.hpp"
#include "anm-lib.hpp"


using namespace std;
using namespace loos;
namespace po = boost::program_options;
using namespace ENM;



// Globals...  Icky poo!

string selection;
string model_name;
string prefix;

int verbosity;
bool debug;

string spring_desc;

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
    "The spring constant used is controlled by the --spring option.\n"
    "If only the name for the spring function is given, then the default\n"
    "parameters are used.  Alternatively, the name may include a\n"
    "comma-separated list of parameters to be passed to the spring\n"
    "function, i.e. --spring=distance,15.0\n\n"
    "Available spring functions:\n";

  vector<string> names = springNames();
  for (vector<string>::iterator i = names.begin(); i != names.end(); ++i)
    cout << "\t" << *i << endl;

  cout <<
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
      ("spring,S", po::value<string>(&spring_desc)->default_value("distance"),"Spring function to use")
      ("fullhelp", "More detailed help");

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
  SpringFunction* spring = 0;
  spring = springFactory(spring_desc);

  SuperBlock* blocker = new SuperBlock(spring, subset);

  ANM anm(blocker);
  anm.debugging(debug);
  anm.prefix(prefix);
  anm.meta(header);
  anm.verbosity(verbosity);

  anm.solve();


  // Write out the LSVs (or eigenvectors)
  writeAsciiMatrix(prefix + "_U.asc", anm.eigenvectors(), header, false);
  writeAsciiMatrix(prefix + "_s.asc", anm.eigenvalues(), header, false);

  writeAsciiMatrix(prefix + "_Hi.asc", anm.inverseHessian(), header, false);

  delete blocker;
  delete spring;
}
