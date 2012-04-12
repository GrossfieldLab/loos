/*
  Matrix cosine content.  Calculates the cosine content 
  of a matrix.  Expected input is the RSVs of a PCA...

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
string rsv;


string fullHelpMessage() {

  string s = 
    "\n"
    "SYNOPSIS\n"
    "\n"
    "Calculate the cosine constent of a right singular vector matrix\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "This tool calculates the cosine content of a matrix.\n"
    "It is intended to be used on the right sinular vectors\n"
    "from an SVD.  These are projections onto to principal\n"
    "components of the simulation.\n"
    "\n"
    "See: Hess, B. \"Convergence of sampling in protein \n"
    "      simulations.\" Phys Rev E (2002) 65(3):031910\n"
    //Matrix cosine content.  Calculates the cosine content 
    //of a matrix.  Expected input is the RSVs of a PCA...
    "\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "rsv-coscon pca_V.asc\n"
    "\tCompute the cos content of the first 10 (default)\n"
    "\tright singular vectors from a simulation PCA.  If\n"
    "\tthe PCA was computed with the LOOS SVD tool, the \n"
    "\tRSVs are stored in _V.asc\n"
    "\n"
    "rsv-coscon --modes=5 pca_V.asc\n"
    "\tCompute the cos content of the first 5 RSVs only.\n"
    "\n"
    "SEE ALSO\n"
    "Packages/Convergence/coscon - \n"
    "\tCompute the cosine content of a matrix.  This tool\n"
    "\tperforms a similar analysis, but it uses a block\n"
    "\taveraging approach where the cosine content is\n"
    "\tcalculated for increasingly long trajectory blocks\n"
    "\n"
    "Packages/Convergence/qcoscon - \n"
    "\tPerform a quick cos content analysis on a simulation.\n"
    "\tSimilar to coscon, but only performs the analysis on\n"
    "\tthe full length simulation.\n"
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


int main(int argc, char *argv[]) {

  string hdr = invocationHeader(argc, argv);
  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  //  opts::BasicOptions* bopts = new opts::BasicOptions;
  ToolOptions* topts = new ToolOptions;
  opts::RequiredArguments* ropts = new opts::RequiredArguments("rsv", "right-singular-vectors");

  opts::AggregateOptions options;
  options.add(bopts).add(topts).add(ropts);
  if (!options.parse(argc, argv))
    exit(-1);

  cout << "# " << hdr << endl;
  cout << "# " << vectorAsStringWithCommas<string>(options.print()) << endl;

  RealMatrix V;
  readAsciiMatrix(ropts->value("rsv"), V);

  cout << "# n\tcoscon\n";
  for (uint i=0; i<nmodes; ++i)
    cout << i << "\t" << cosineContent(V, i) << endl;

}
