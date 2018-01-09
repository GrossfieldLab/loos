/*
  distance_prob.cpp

  (c) 2018 Alan Grossfield
           Department of Biochemistry
           University of Rochster School of Medicine and Dentistry

  Compute atom-atom distance probability function


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
// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  double hist_min, hist_max;
  int num_bins;
  string prefix;
  bool use_electrons;
  bool write_per_frame;

  // Change these options to reflect what your tool needs
  void addGeneric(po::options_description& o) {
    o.add_options()
    ("hist_min", po::value<double>(&hist_min)->default_value(0.0), "Histogram minimum")
    ("hist_max", po::value<double>(&hist_max)->default_value(50.0), "Histogram maximum")
    ("num_bins", po::value<int>(&num_bins)->default_value(100), "Number of bins")
    ("prefix", po::value<string>(&prefix)->default_value(string("./foo_")), "Output file prefix")
    ("electrons", "Weight atoms by electrons")
    ("per-frame", "Write a distribution for each frame")
    ;
  }

    bool postConditions(po::variables_map& vm)
    {
      if (vm.count("electrons")) {
        use_electrons = true;
      }
      else {
        use_electrons = false;
      }
      if (vm.count("per-frame")) {
        write_per_frame = true;
      }
      else {
        write_per_frame = false;
      }
      return true;
    }


  // The print() function returns a string that describes what all the
  // options are set to (for logging purposes)
  string print() const {
    ostringstream oss;
    oss << boost::format("hist_min=%f, hist_max=%f, num_bins=%d, prefix=%s")
            % hist_min
            % hist_max
            % num_bins
            % prefix.c_str();
    return(oss.str());
  }

};
// @endcond
// ----------------------------------------------------------------

string fullHelp(void) {
  string s =
  "\n"
  "SYNOPSIS\n"
  "\n"
  "Compute electron-weighted atom-atom distance distribution function\n"
  "\n"
  "DESCRIPTION\n"
  "\n"
  "This tool is designed to produce a pair-distribution function \n"
  "comparable to what you'd get from an X-ray scattering experiment.\n"
  "Given a selection, for each frame it computes the pair distance \n"
  "distribution function, either weighting each pair equally \n"
  "or by the product of their number of electrons.\n"
  "\n"
  "WARNING: this means you need charge and mass information (to deduce \n"
  "the atomic number).  If you use something other than a PSF to define \n"
  "the system, this information won't be available, and you'll \n"
  "probably get the unweighted distance distribution function instead.\n"
  "\n"
  "For the moment, it is hardwired to write out a distribution for each \n"
  "frame.\n"
  "\n";
  return(s);
}



int main(int argc, char *argv[]) {

  // Store the invocation information for logging later
  string header = invocationHeader(argc, argv);

  // Build up the command-line options for this tool
  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelp());
  opts::BasicSelection* sopts = new opts::BasicSelection;
  opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
  ToolOptions* topts = new ToolOptions;

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
  vector<double> weighting;
  weighting.assign(subset.size(), 0.0);
  if (topts->use_electrons) {
    subset.deduceAtomicNumberFromMass();
    for (uint i=0; i<subset.size(); i++) {
      weighting[i] = subset[i]->atomic_number() - subset[i]->charge();
    }
  }
  else {
    weighting.assign(subset.size(), 1.0);
  }

  vector<double> total_histogram;
  total_histogram.assign(topts->num_bins, 0.0);
  uint frames_accumulated = 0;

  // Now iterate over all frames in the trajectory (excluding the skip
  // region)
  while (traj->readFrame()) {
    // Set up the histogram
    vector<double> histogram;
    histogram.assign(topts->num_bins, 0.0);
    double normalization = 0.0;
    int excluded = 0;

    // Update the coordinates ONLY for the subset of atoms we're
    // interested in...
    traj->updateGroupCoords(subset);
    GCoord box = model.periodicBox();

    for (uint i=0; i<subset.size()-1; i++) {
      for (uint j=i+1; j<subset.size(); j++) {
        double distance = subset[i]->coords().distance(subset[j]->coords(),
                                                         box);
        if (distance >= topts->hist_max) {
          excluded++;
        }
        else {
         int bin = static_cast<int>((distance - topts->hist_min)/bin_width);
         double e2 = weighting[i] * weighting[j];
         histogram.at(bin) += e2;
         normalization += e2;
        }
      }
    }

    if (excluded) {
      cerr << "Frame: " << traj->currentFrame()
           << " excluded " << excluded
           << " distances." << endl;
    }

    // Output the histogram for the frame
    if (topts->write_per_frame) {
      string filename = topts->prefix + to_string(traj->currentFrame())
                                      + string(".dat");
      ofstream outfile(filename.c_str());
      if (outfile.fail()) {
        cerr << "Error opening file " << filename << endl;
        exit(-1);
      }
      outfile << "# Distance Probability" << endl;
      for (uint i=0; i<histogram.size(); i++) {
        histogram[i] /= normalization;
        double d = topts->hist_min + (i+0.5)*bin_width;
        outfile << d << "\t" << histogram[i] << endl;
      }
    }

    // Accumulate the total histogram
    for (uint i=0; i<histogram.size(); i++) {
      total_histogram[i] += histogram[i];
      frames_accumulated++;
    }


  }

  cout << "# Distance Probability" << endl;
  for (uint i=0; i<total_histogram.size(); i++) {
    total_histogram[i] /= frames_accumulated;
    double d = topts->hist_min + (i+0.5)*bin_width;
    cout << d << "\t" << total_histogram[i] << endl;
  }
}
