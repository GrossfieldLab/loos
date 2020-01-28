/*
  Compute the distribution or timeseries of radii of gyration for a selection of
  atoms.


  Alan Grossfield
  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2020, Louis Smith
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

#include <iostream>
#include <loos.hpp>

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

const string fullHelpMessage = "XXX";

class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() {}
  // clang-format off
  void addGeneric(po::options_description& o) { 
    o.add_options()
      ("timeseries,t", po::value<string>(&timeseries)->default_value(""), "Write frame-by-frame timeseries to file name provided. If none provided, not written.")
      ("num-bins,n", po::value<int>(&num_bins)->default_value(50), "Number of bins to use for histogramming.")
      ("min-bin,m", po::value<double>(&min_bin)->default_value(0), "Minimum value for the histogram bins.")
      ("max-bin,M", po::value<double>(&max_bin)->default_value(50), "Maximum value for the histogram bins")
      ("by-molecule", po::value<bool>(&by_molecule)->default_value(false), "Split provided selection by connectivity within that selection." )
    ;
  }
  // clang-format on
  string print() const {
    ostringstream oss;
    oss << boost::format("timeseries=%s,min_bin=%s,max_bin=%s,by_molecule=%b,"
                         "num_bins=%d") %
               timeseries % min_bin % max_bin % by_molecule % num_bins;
    return (oss.str());
  }
  string timeseries;
  double min_bin;
  double max_bin;
  bool by_molecule;
  int num_bins;
};

inline void histogram_rgyr(vector<greal> &hist, const greal rgyr,
                           const greal min_bin, const greal max_bin,
                           const greal bin_width, int &count, const int frame,
                           ofstream &outfile) {
  if ((rgyr >= min_bin) && (rgyr < max_bin)) {
    hist[int((rgyr - min_bin) / bin_width)]++;
    count++;
  }
}

inline void ts_hist_rgyr(vector<greal> &hist, const greal rgyr,
                         const greal min_bin, const greal max_bin,
                         const greal bin_width, int &count, const int frame,
                         ofstream &outfile) {
  if ((rgyr >= min_bin) && (rgyr < max_bin)) {
    hist[int((rgyr - min_bin) / bin_width)]++;
    count++;
  }
  outfile << frame << "\t" << rgyr << "\n";
}

int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);

  opts::BasicOptions *bopts = new opts::BasicOptions(fullHelpMessage);
  opts::BasicSelection *sopts = new opts::BasicSelection("all");
  opts::MultiTrajOptions *mtopts = new opts::MultiTrajOptions;
  ToolOptions *topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(mtopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  // Write histogram results to stdout. Write timeseries, if asked for, to file
  cout << "# " << header << "\n";
  ofstream tsf(topts->timeseries);
  // make a function pointer with a signature matching the
  void (*frameOperator)(vector<greal> & hist, const greal rgyr,
                        const greal min_bin, const greal max_bin,
                        const greal bin_width, int &count, const int frame,
                        ofstream &outfile);
  // pick which operation to perform per frame using function pointer
  if (topts->timeseries.empty())
    frameOperator = histogram_rgyr;
  else {
    tsf << "# " << header << "\n"
        << "# Rgyr\tProb\tCum\n";
    frameOperator = ts_hist_rgyr;
  }

  // establish system, and molecular subsystems
  vector<AtomicGroup> molecules;
  if (topts->by_molecule)
    molecules = mtopts->model.splitByMolecule(sopts->selection);
  else
    molecules.push_back(selectAtoms(mtopts->model, sopts->selection));

  // prepare for trajectory loop
  // counter for number of molecules in histogram bounds
  int count = 0;
  const int num_bins = topts->num_bins;
  const greal min_bin = topts->min_bin;
  const greal max_bin = topts->max_bin;
  const greal bin_width = (max_bin - min_bin) / num_bins;

  // define and zero histogram
  vector<greal> hist(num_bins, 0.0);
  // place to put the radius of gyration of each molecule
  greal rgyr;
  while (mtopts->trajectory->readFrame()) {
    mtopts->trajectory->updateGroupCoords(mtopts->model);
    for (AtomicGroup mol : molecules) {
      rgyr = mol.radiusOfGyration();
      // call function pointer, passing in all the state from the loop needed
      // for either histogramming or timeseries writing
      (*frameOperator)(hist, rgyr, min_bin, max_bin, bin_width, count,
                       mtopts->trajectory->currentFrame(), tsf);
    }
  }
  // Write the histogram to stdout
  cout << "# Rgyr\tProb\tCum" << endl;
  greal cum = 0.0;
  for (int i = 0; i < num_bins; i++) {
    greal d = bin_width * (i + 0.5) + min_bin;

    greal prob = hist[i] / count;
    cum += prob;
    cout << d << "\t" << prob << "\t" << cum << endl;
  }
  tsf.close();
}