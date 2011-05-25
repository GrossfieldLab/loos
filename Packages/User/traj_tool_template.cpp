/*
  traj_tool_template.cpp

  (c) 2011 Tod D. Romo, Grossfield Lab
           Department of Biochemistry
           University of Rochster School of Medicine and Dentistry


  C++ template for writing a tool that works with a trajectory
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2011 Tod D. Romo
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



// ----------------------------------------------------------------
// The following code is for implementing tool-specific
// options if the tool requires them.  If not, this section can be
// deleted.

double option1;
int option2;


// The following conditional prevents this class from appearing in the
// Doxygen documentation for LOOS:
//
// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:

  // Change these options to reflect what your tool needs
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("option1", po::value<double>(&option1)->default_value(0.0), "Tool Option #1")
      ("option2", po::value<int>(&option2)->default_value(42), "Tool option #2");
  }

  // The print() function returns a string that describes what all the
  // options are set to (for logging purposes)
  string print() const {
    ostringstream oss;
    oss << boost::format("option1=%f, option2=%d") % option1 % option2;
    return(oss.str());
  }

};
// @endcond
// ----------------------------------------------------------------


void calculate(const AtomicGroup& structure) {
  // Do something here with an atom...
}




int main(int argc, char *argv[]) {
  
  // Store the invocation information for logging later
  string header = invocationHeader(argc, argv);
  
  // Build up the command-line options for this tool by instantiating
  // the appropriate OptionsPackage objects...

  // Basic options should be used by all tools.  It provides help,
  // verbosity, and the ability to read options from a config file
  opts::BasicOptions* bopts = new opts::BasicOptions;

  // This tool can operate on a subset of atoms.  The BasicSelection
  // object provides the "--selection" option.
  opts::BasicSelection* sopts = new opts::BasicSelection;

  // The BasicTrajectory object handles specifying a trajectory as
  // well as a "--skip" option that lets the tool skip the first
  // number of frames (i.e. equilibration).  It creates a pTraj object
  // that is already primed for reading...
  opts::BasicTrajectory* tropts = new opts::BasicTrajectory;

  // Tool-specific options can be included here...
  ToolOptions* topts = new ToolOptions;

  // All of the OptionsPackages are combined via the AggregateOptions
  // object.  First instantiate it, then add the desired
  // OptionsPackage objects.  The order is important.  We recommend
  // you progress from general (Basic and Selection) to more specific
  // (model) and finally the tool options.
  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(tropts).add(topts);

  // Parse the command-line.  If an error occurred, help will already
  // be displayed and it will return a FALSE value.
  if (!options.parse(argc, argv))
    exit(-1);

  // Pull the model from the options object (it will include coordinates)
  AtomicGroup model = tropts->model;
  
  // Pull out the trajectory...
  pTraj traj = tropts->trajectory;

  // Select the desired atoms to operate over...
  AtomicGroup subset = selectAtoms(model, sopts->selection);

  // Now iterate over all frames in the trajectory (excluding the skip
  // region)
  while (traj->readFrame()) {

    // Update the coordinates ONLY for the subset of atoms we're
    // interested in...
    traj->updateGroupCoords(subset);

    // Now calculate something with the AtomicGroup
    calculate(subset);
  }

}
