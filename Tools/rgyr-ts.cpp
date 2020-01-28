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

#include <loos.hpp>
#include <iostream>

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
      ("num_bins,n", po::value<int>(&num_bins)->default_value(50), "Number of bins to use for histogramming.")
      ("bin-min,m", po::value<double>(&bin_min)->default_value(0), "Minimum value for the histogram bins.")
      ("bin-max,M", po::value<double>(&bin_max)->default_value(50), "Maximum value for the histogram bins")
      ("by-molecule", po::value<bool>(&by_molecule)->default_value(false), "Split provided selection by connectivity within that selection." )
    ;
  }
// clang-format on
  string print() const {
    ostringstream oss;
    oss << boost::format("timeseries=%s") % timeseries;
    return(oss.str());
  }
  string timeseries;
  double bin_min;
  double bin_max;
  bool by_molecule;
  int num_bins;
};

void histogram_rgyr(vector<greal>& hist, greal rgyr, greal hist_min, greal hist_max, greal bin_width, int count, int frame, ofstream& outfile) {
  if ((rgyr >= hist_min) && (rgyr < hist_max)){
    hist[int((rgyr-hist_min)/bin_width)]++;
    count++;
  }
}

void ts_hist_rgyr(vector<greal>& hist, greal rgyr, greal hist_min, greal hist_max, greal bin_width, int count, int frame,  ofstream& outfile){
  if ((rgyr >= hist_min) && (rgyr < hist_max)){
    hist[int((rgyr-hist_min)/bin_width)]++;
    count++;
  }
  outfile << frame << "\t" << rgyr << "\n";
}

int main (int argc, char * argv[]) {
  string header = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage);
  opts::BasicSelection* sopts = new opts::BasicSelection("all");
  opts::MultiTrajOptions* mtopts = new opts::MultiTrajOptions;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(mtopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  // Write histogram results to stdout. Write timeseries, if asked for, to file
  cout << "# " << header << "\n";
  ofstream tsf(topts->timeseries);
  void (*frameOperator)(vector<greal>& hist, greal rgyr, greal hist_min, greal hist_max, greal bin_width, int count, int frame,  ofstream& outfile);
  // pick which operation to perform per frame using function pointer
  if(topts->timeseries.empty())
    frameOperator = histogram_rgyr;
  else
    frameOperator = ts_hist_rgyr;    
  
  int 
}