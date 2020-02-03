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

const string fullHelpMessage =
    // clang-format off
    "SYNOPSIS \n"
    " \n"
    "Read a set of trajectories and return a histogram of the radius of \n"
    "gyration of a selection. \n"
    " \n"
    "DESCRIPTION \n"
    " \n"
    "This tool computes the radius of gyration for a selection of atoms for each \n"
    "frame of a provided multi-trajectory. By default this selection is treated as \n"
    "one group, but if --by-molecule is thrown then the selection will be split by \n"
    "connectivity. The skip, stride, and range options operate on the multi-\n"
    "trajectory as they would for other LOOS tools. The time-series option specifies\n"
    "a file to write a timeseries of the radius of gyrations to. If none is \n"
    "specified then no time-series is written. The num-bins, min-bin, and max-bin \n"
    "options determine the extent and bin-width of the histogram, which is written \n"
    "to stdout. The histogram contains both the probability per bin, and the \n"
    "cumulative probability, in three tab-delimited columns with the leftmost being\n"
    "the bin-center. \n"
    " \n"
    "EXAMPLE \n"
    " \n"
    "rad-gyr -k 100 -n 20 -m 5 -M 25 --by-molecule model.pdb traj1.dcd traj2.dcd \\\n"
    "traj3.dcd \n"
    " \n"
    "This will concatenate traj1, traj2, and traj3 into one virtual trajectory, skip\n"
    "the first 100 frames of each, and then histogram the radius of gyration of all\n"
    "the atoms in each of the molecules in the model provided with the \n"
    "trajectories, computing their radii of gyration and summing over them. Note \n"
    "that this would be nonsensical if the molecules in the model were not multiple \n"
    "copies of the same molecule, since all such Rgyr will be collated into one \n"
    "histogram at the end. Using selection to pick a subsystem of interest makes \n"
    "more sense in the context of most analyses of a solute.\n";
//clang-format on
class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() {}
  // clang-format off
  void addGeneric(po::options_description& o) { 
    o.add_options()
      ("timeseries,t", po::value<string>(&timeseries)->default_value(""), 
       "Write Rgyr per-frame to file name provided.")
      ("num-bins,n", po::value<int>(&num_bins)->default_value(50), 
       "Number of bins to use for histogramming.")
      ("min-dist,m", po::value<double>(&min_dist)->default_value(0), 
       "Minimum value for the histogram bins.")
      ("max-dist,M", po::value<double>(&max_dist)->default_value(50), 
       "Maximum value for the histogram bins.")
      ("by-molecule", po::bool_switch(&by_molecule)->default_value(false), 
       "Split 'selection' by connectivity of 'model'.")
      ("by-fragment", po::bool_switch(&by_molecule)->default_value(false), 
       "Split 'selection' by its own connectivity.")
       
    ;
  }
  // clang-format on
  string print() const {
    ostringstream oss;
    oss << boost::format("timeseries=%s,min_dist=%s,max_dist=%s,by_molecule=%b,"
                         "by_fragment=%b,num_bins=%d") %
               timeseries % min_dist % max_dist % by_molecule % by_fragment %
               num_bins;
    return (oss.str());
  }

  bool postConditions(po::variables_map &map) {
    if(by_molecule && by_fragment){
      cerr << "ERROR: --by-molecule and --by-fragment flags are mutually exclusive.\n";
      return(false);
    } else
      return(true);
    
  }
  string timeseries;
  double min_dist;
  double max_dist;
  bool by_molecule;
  bool by_fragment;
  int num_bins;
};

// Do the work of recording, either in the case where a TS is requested,
inline void histogram_rgyr(vector<greal> &hist, const greal rgyr,
                           const greal min_dist, const greal max_dist,
                           const greal bin_width, int &count,
                           ofstream &outfile) {
  if ((rgyr >= min_dist) && (rgyr < max_dist)) {
    hist[int((rgyr - min_dist) / bin_width)]++;
    count++;
  }
}
inline void histogram_molecules_rgyr(vector<greal> &hist,
                                     vector<AtomicGroup> &molecules,
                                     const greal min_dist, const greal max_dist,
                                     const greal bin_width, int &count,
                                     const int frame, ofstream &outfile) {
  for (auto molecule : molecules) {
    histogram_rgyr(hist, molecule.radiusOfGyration(), min_dist, max_dist,
                   bin_width, count, outfile);
  }
}
// or in the case when not.
inline void ts_hist_rgyr(vector<greal> &hist, vector<AtomicGroup> &molecules,
                         const greal min_dist, const greal max_dist,
                         const greal bin_width, int &count, const int frame,
                         ofstream &outfile) {
  greal rgyr;
  outfile << frame;
  for (auto molecule : molecules) {
    rgyr = molecule.radiusOfGyration();
    outfile << "\t" << rgyr;
    histogram_rgyr(hist, rgyr, min_dist, max_dist, bin_width, count, outfile);
  }
  outfile << "\n";
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
  void (*frameOperator)(vector<greal> & hist, vector<AtomicGroup> & molecules,
                        const greal min_dist, const greal max_dist,
                        const greal bin_width, int &count, const int frame,
                        ofstream &outfile);

  // establish system, and molecular subsystems
  vector<AtomicGroup> molecules;
  if (topts->by_molecule)
    molecules = mtopts->model.splitByMolecule(sopts->selection);
  else if(topts->by_fragment)
    molecules = selectAtoms(mtopts->model, sopts->selection).splitByMolecule();
  else
    molecules.push_back(selectAtoms(mtopts->model, sopts->selection));

  // pick which operation to perform per frame using function pointer
  if (topts->timeseries.empty())
    frameOperator = histogram_molecules_rgyr;
  else {
    tsf << "# " << header << "\n"
        << "# frame";
    // Label each column with the starting and ending index for the 'molecule'
    // as a stand in for the case where molecule's chainids are inaccurate.
    for (auto molecule : molecules) {
      int first = molecule[0]->index();
      int last = (*(molecule.end() - 1))->index();
      tsf << "\tatoms" << first << "-" << last;
    }
    tsf << "\n";
    frameOperator = ts_hist_rgyr;
  }
  // prepare for trajectory loop
  // counter for number of molecules in histogram bounds
  int count = 0;
  const int num_bins = topts->num_bins;
  const greal min_dist = topts->min_dist;
  const greal max_dist = topts->max_dist;
  const greal bin_width = (max_dist - min_dist) / num_bins;

  // define and zero histogram
  vector<greal> hist(num_bins, 0.0);
  while (mtopts->trajectory->readFrame()) {
    mtopts->trajectory->updateGroupCoords(mtopts->model);
    // call fxn pointer to accumulate histogram and optionally write timeseries
    (*frameOperator)(hist, molecules, min_dist, max_dist, bin_width, count,
                     mtopts->trajectory->currentFrame(), tsf);
  }
  // Write the histogram to stdout
  cout << "# Rgyr\tProb\tCum" << endl;
  greal cum = 0.0;
  for (int i = 0; i < num_bins; i++) {
    greal d = bin_width * (i + 0.5) + min_dist;

    greal prob = hist[i] / count;
    cum += prob;
    cout << d << "\t" << prob << "\t" << cum << endl;
  }
  tsf.close();
}