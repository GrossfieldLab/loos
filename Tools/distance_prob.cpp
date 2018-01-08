/*
  traj_calc.cpp

  (c) 2011 Tod D. Romo, Grossfield Lab
           Department of Biochemistry
           University of Rochster School of Medicine and Dentistry


  C++ template for writing a tool that performs a calculation on a trajectory
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
// Tool-specific options



// The following conditional prevents this class from appearing in the
// Doxygen documentation for LOOS:
//
// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  double hist_min, hist_max;
  int num_bins;
  string prefix;

  // Change these options to reflect what your tool needs
  void addGeneric(po::options_description& o) {
    o.add_options()
    ("hist_min", po::value<double>(&hist_min)->default_value(0.0), "Histogram minimum")
    ("hist_max", po::value<double>(&hist_min)->default_value(50.0), "Histogram maximum")
    ("num_bins", po::value<int>(&num_bins)->default_value(100), "Number of bins")
    ("prefix", po::value<string>(&prefix)->default_value(string("./foo_")), "Output file prefix")
    ;
  }

  /*
  // The print() function returns a string that describes what all the
  // options are set to (for logging purposes)
  string print() const {
    ostringstream oss;
    oss << boost::format("option1=%f, option2=%d") % option1 % option2;
    return(oss.str());
  }
  */

};
// @endcond
// ----------------------------------------------------------------




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

  // ***EDIT***
  // Tool-specific options can be included here...
  ToolOptions* topts = new ToolOptions;

  // ***EDIT***
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

  AtomicGroup model = tropts->model;
  pTraj traj = tropts->trajectory;

  // Select the desired atoms to operate over...
  AtomicGroup subset = selectAtoms(model, sopts->selection);

  // For convenience, figure out histogram bin width
  double bin_width = (topts->hist_max - topts->hist_min)/topts->num_bins;

  // Compute electrons for each atom
  subset.deduceAtomicNumberFromMass();
  vector<double> electrons(subset.size());
  for (uint i=0; i<subset.size(); i++) {
    electrons[i] = subset[i]->atomic_number() - subset[i]->charge();
  }


  // Now iterate over all frames in the trajectory (excluding the skip
  // region)
  while (traj->readFrame()) {
    // Set up the histogram
    vector<double> histogram(topts->num_bins);
    double normalization = 0.0;

    // Update the coordinates ONLY for the subset of atoms we're
    // interested in...
    traj->updateGroupCoords(subset);
    GCoord box = model.periodicBox();

    for (uint i=0; i<subset.size()-1; i++) {
       for (uint j=i+1; j<subset.size(); j++) {
         double distance = subset[i]->coords().distance(subset[j]->coords(),
                                                        box);
         int bin = static_cast<int>((distance - topts->hist_min)/bin_width);
         double e2 = electrons[i] * electrons[j];
         histogram[bin] += e2;
         normalization += e2;
       }
    }

    // Output the histogram for the frame

    string filename = topts->prefix + to_string(traj->currentFrame());
    ofstream outfile(filename.c_str());
    outfile << "# Distance Probability" << endl;
    for (uint i=0; i<histogram.size(); i++) {
      histogram[i] /= normalization;
      double d = topts->hist_min + (i+0.5)*bin_width;
      outfile << d << "\t" << histogram[i] << endl;
    }

  }

}
