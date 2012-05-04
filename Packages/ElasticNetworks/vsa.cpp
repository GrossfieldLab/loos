/*
  vsa

  Usage:
    vsa [options] subset environment model output_prefix

  See:
    Woodcock et al, J Chem Phys (2008) 129:214109
    Haffner & Zheng, J Chem Phys (2009) 130:194111

*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009 Tod D. Romo
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
#include "vsa-lib.hpp"


using namespace std;
using namespace loos;
using namespace ENM;


namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;


// Globals...
double normalization = 1.0;
double threshold = 1e-10;

string hdr;
string subsystem_selection, environment_selection, model_name, prefix;
int verbosity = 0;
bool debug = false;
bool occupancies_are_masses;
string psf_file;

string spring_desc;
bool nomass;


string fullHelpMessage() {

  string s =
    "\n"
    "SYNOPSIS\n"
    "\n"
    "Compute the vibrational subsystem analysis version of the\n"
    "anisotropic network model.  (See Woodcock et al.)\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "Computes the VSA network model given a subsystem and an\n"
    "environment selection.  The output consists of several different\n"
    "ASCII formatted matrices (that can be read by Matlab/Octave) and\n"
    "depends on whether or not masses are included in the\n"
    "calculation.  If debugging is turned on (--debug), then the\n"
    "intermediate matrices are written out:\n"
    "\tfoo_H.asc    - Composite Hessian\n"
    "\tfoo_Hss.asc  - Subsystem Hessian\n"
    "\tfoo_Hee.asc  - Environment Hessian\n"
    "\tfoo_Hse.asc  - Subsystem-Environment Hessian\n"
    "\tfoo_Heei.asc - Inverted Environment Hessian\n"
    "\tfoo_Hssp.asc - Effective Subsystem Hessian\n"
    "\tfoo_Ms.asc   - Subsystem mass (optional)\n"
    "\tfoo_Me.asc   - Environment mass (optional)\n"
    "\tfoo_Msp.asc  - Effective subsystem mass (optional)\n"
    "\tfoo_R.asc    - Cholesky decomposition of Msp (optional)\n"
    "\n\n"
    "* Unit Subsystem Mass, Zero Environment Mass *\n\n"
    "Here, the effective subsystem Hessian is created and a Singular\n"
    "Value Decomposition used to solve the eigenproblem:\n"
    "\tfoo_U.asc = Subsystem eigenvectors\n"
    "\tfoo_s.asc = Subsystem eigenvalues\n"
    "\n\n"
    "* Subsystem and Environment Mass *\n\n"
    "The generalized eigenvalue problem is solved creating the\n"
    "following matrices:\n"
    "\tfoo_s.asc = Subsystem eigenvalues (mass-weighted)\n"
    "\tfoo_U.asc = Subsystem eigenvectors (mass-weighted)\n"
    "\n\n"
    "* Spring Constant Control *\n\n"
    "The spring constant used is controlled by the --spring option.\n"
    "If only the name for the spring function is given, then the default\n"
    "parameters are used.  Alternatively, the name may include a\n"
    "comma-separated list of parameters to be passed to the spring\n"
    "function, i.e. --spring=distance,15.0\n\n"
    "Available spring functions:\n";

  vector<string> names = springNames();
  for (vector<string>::iterator i = names.begin(); i != names.end(); ++i)
    s = s + "\t" + *i + "\n";

  s +=
    "\n\n"
    "* Mass Control *\n\n"
    "VSA, by default, assumes that masses will be present.  These can\n"
    "come from one of two sources.  If \"--psf foo.psf\" is given,\n"
    "then the masses will be assigned using the \"foo.psf\" file.  This\n"
    "assumes that the atoms are in the same order between the PSF file\n"
    "and the structure file given on the command line.  Alternatively,\n"
    "the occupancy field of the PDB can be used with the\n"
    "\"--occupancies 1\" option.  See the psf-masses tool for one way to\n"
    "copy masses into a PDB's occupancies.\n"
    "\n"
    "To disable masses (i.e. use unit masses for the subsystem and\n"
    "zero masses for the environment), use the \"--nomass 1\" option.\n"
    "\n\n"
    "EXAMPLES \n\n"
    "\n"
    "vsa --occupancies 1 foo.pdb 'segid == \"TRAN\" && name == \"CA\"'\\\n"
    "  'segid != \"TRAN\" && name == \"CA\"' foo_vsa\n"
    "\tCompute the VSA for a transmembrane region based on segid with the\n"
    "\tmasses stored in the occupancy field of the PDB.  Here the  enviroment\n"
    "\tcontains all other CA's in the system.\n"
    "\n"
    "vsa --psf foo.psf foo.pdb \"`cat selection` && name == 'CA'\" \\\n"
    "   \"not (`cat selection`) && name == 'CA'\" foo_vsa\n"
    "\tCompute the VSA for a transmembrane region where the selection\n"
    "\tis stored in a file and masses are taken from a PSF file.\n"
    "\n"
    "vsa --nomass 1 foo.pdb 'name == \"CA\"' 'name =~ \"^(C|O|N)$\"' foo_vsa\n"
    "\tCompute the mass-less VSA with CAs as the subsystem and all other\n"
    "\tbackbone atoms as the environment.\n"
    "\n"
    "vsa --nomass 1 --spring hca foo.pdb 'name == \"CA\"' 'name =~ \"^(C|O|N)$\"' foo_vsa\n"
    "\tThe same example as above, but using the HCA spring constants.\n";

  return(s);
    
}




// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("psf", po::value<string>(&psf_file), "Take masses from the specified PSF file")
      ("debug", po::value<bool>(&debug)->default_value(false), "Turn on debugging (output intermediate matrices)")
      ("occupancies", po::value<bool>(&occupancies_are_masses)->default_value(false), "Atom masses are stored in the PDB occupancy field")
      ("nomass", po::value<bool>(&nomass)->default_value(false), "Disable mass as part of the VSA solution")
      ("spring", po::value<string>(&spring_desc)->default_value("distance"), "Spring method and arguments");
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("psf='%s', debug=%d, occupancies=%d, nomass=%d, spring='%s'")
      % psf_file
      % debug
      % occupancies_are_masses
      % nomass
      % spring_desc;
    return(oss.str());
  }


};
// @endcond




int main(int argc, char *argv[]) {
  hdr = invocationHeader(argc, argv);

  // Build up options
  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::ModelWithCoords* mopts = new opts::ModelWithCoords;
  ToolOptions* topts = new ToolOptions;
  opts::RequiredArguments* ropts = new opts::RequiredArguments;
  ropts->addArgument("subsystem", "subsystem-selection");
  ropts->addArgument("environment", "environment-selection");
  ropts->addArgument("prefix", "output-prefix");

  opts::AggregateOptions options;
  options.add(bopts).add(mopts).add(topts).add(ropts);
  if (!options.parse(argc, argv))
    exit(-1);

  // Extract values
  AtomicGroup model = mopts->model;
  verbosity = bopts->verbosity;
  subsystem_selection = ropts->value("subsystem");
  environment_selection = ropts->value("environment");
  prefix = ropts->value("prefix");

  // Ugly way of handling multiple methods for getting masses into the equation...
  if (verbosity > 0)
    cerr << "Assigning masses...\n";

  // Get masess from somewhere, or rely on the defaults provided by
  // Atom (i.e. should be 1)
  if (! psf_file.empty())
    massFromPSF(model, psf_file);
  else if (occupancies_are_masses)
    massFromOccupancy(model);
  else if (!nomass)
    cerr << "WARNING- using default masses\n";


  // Partition the model for building the composite Hessian
  AtomicGroup subsystem = selectAtoms(model, subsystem_selection);
  AtomicGroup environment = selectAtoms(model, environment_selection);
  AtomicGroup composite = subsystem + environment;

  if (verbosity > 1) {
    cerr << "Subsystem size is " << subsystem.size() << endl;
    cerr << "Environment size is " << environment.size() << endl;
  }

  // Determine which kind of scaling to apply to the Hessian...
  SpringFunction* spring = 0;
  spring = springFactory(spring_desc);

  SuperBlock* blocker = new SuperBlock(spring, composite);

  VSA vsa(blocker, subsystem.size());
  vsa.prefix(prefix);
  vsa.meta(hdr);
  vsa.debugging(debug);
  vsa.verbosity(verbosity);

  if (!nomass) {
    DoubleMatrix M = getMasses(composite);
    vsa.setMasses(M);
  }

  vsa.solve();
  
  writeAsciiMatrix(prefix + "_U.asc", vsa.eigenvectors(), hdr, false);
  writeAsciiMatrix(prefix + "_s.asc", vsa.eigenvalues(), hdr, false);

  // Be good...
  delete spring;
  delete blocker;

}
