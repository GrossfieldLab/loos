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
  opts::BasicOptions* bopts = new opts::BasicOptions;
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
