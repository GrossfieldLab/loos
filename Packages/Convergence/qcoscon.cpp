/*
  
  "Quick" Cosine content.  Calculates the cosine content for the
  entire trajectory for the first N modes.

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


uint nmodes;

string fullHelpMessage() {

  string s = 
    "\n"
    "SYNOPSIS\n"
    "\n"
    "Calculate the cosine constent of a whole simulation\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "Quick version of the cosine content calculation.\n"
    "This tool performs the same calculation as coscon, \n"
    "but instead of varying the trajectory in a block\n"
    "averaging approach only the full trajectory is used.\n"
    "The results are printed for the first 10 modes.\n"
    "\n"
    //
    "EXAMPLES\n"
    "\n"
    "qcoscon -s 'name==\"CA\"' model.pdb traj.dcd\n"
    "\tCalculate the cos content of the first 10 modes\n"
    "\tof traj.dcd using the PCA of the CA atoms.\n"
    "\n"
    "SEE ALSO\n"
    "Packages/Convergence/coscon - \n"
    "\tCompute the cosine content of a matrix.  This tool\n"
    "\tperforms a similar analysis, but it uses a block\n"
    "\taveraging approach where the cosine content is\n"
    "\tcalculated for increasingly long trajectory blocks\n"
    "\n"
    "Packages/Convergence/rsv-coscon - \n"
    "\tCalculate the cos content of the RSVs from a simulation\n"
    "\tPCA.\n"
    "\n"
    "Tools/svd - \n"
    "\tCompute the principal components via the SVD.\n"
    "\tThis results in several matrix files including\n"
    "\tthe RSVs used as input to the current tool. \n"
    "\tThe file [prefix]_V.asc contains the RSV matrix.\n"
    "\n"
    "\n";

  return(s);
}


//@cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("modes", po::value<uint>(&nmodes)->default_value(10), "Compute cosine content for first N modes");
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("modes=%d") % nmodes;
    return(oss.str());
  }

};


// @endcond TOOLS_INTERNAL

int main(int argc, char *argv[]) {

  string hdr = invocationHeader(argc, argv);
  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
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

  AtomicGroup subset = selectAtoms(model, sopts->selection);
  vector<AtomicGroup> ensemble;
  readTrajectory(ensemble, subset, traj);
 
  // Read in and align...
  boost::tuple<std::vector<XForm>, greal, int> ares = iterativeAlignment(ensemble);
  AtomicGroup avg = averageStructure(ensemble);

  NoAlignPolicy policy(avg, true);
  RealMatrix V = rsv(ensemble, policy);
  cout << "# n\tcoscon\n";
  for (uint i=0; i<nmodes; ++i)
    cout << i << "\t" << cosineContent(V, i) << endl;

}
