/*
  Vertical (along Z) area profile using radius of gyration or max radius
 
  This file is part of LOOS.
  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2016, Tod D. Romo, Grossfield Lab,
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

#include <boost/algorithm/string.hpp>

using namespace std;
using namespace loos;



namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;


// @cond TOOL_INTERNAL
string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "Compute vertical area profile using radius of gyration or max radius\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\n"
    "\n"
    "SEE ALSO\n"
    "\tarea_profile.py\n";
  
  return(msg);
}



class ToolOptions : public opts::OptionsPackage {
public:

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("mode", po::value<string>(&mode_selection)->default_value("rgyr"), "Calculation type (rgyr, maxr)")
      ("tseries", po::value<string>(&timeseries_filename), "Output time series to this file");
  }


  void addHidden(po::options_description& o) {
    o.add_options()
      ("zmin", po::value<double>(&zmin), "Min position along Z")
      ("zmax", po::value<double>(&zmax), "Max position along Z")
      ("nbins", po::value<uint>(&nbins), "Number of bins along Z");
  }


  void addPositional(po::positional_options_description& o) {
    o.add("zmin", 1);
    o.add("zmax", 1);
    o.add("nbins", -1);
  }

  
  bool postConditions(po::variables_map& vm) 
  {

    boost::algorithm::to_lower(mode_selection);
    if (mode_selection == "rgyr")
      rgyr_mode = true;
    else if (mode_selection == "maxr")
      rgyr_mode = false;
    else {
      cerr << "Error- mode must be one of: rgyr, maxr\n";
      return(false);
    }
    
    return(true);
  }


  string help() const {
    return("zmin zmax nbins");
  }
    
  string print() const {
    ostringstream oss;

    oss << boost::format("mode='%s',tseries='%s',zmin=%f,zmax=%f,nbins=%u")
      % mode_selection
      % timeseries_filename
      % zmin
      % zmax
      % nbins;
    return(oss.str());
  }


  string mode_selection, timeseries_filename;
  bool rgyr_mode;
  double zmin, zmax;
  uint nbins;
};


typedef vector<AtomicGroup>   Slices;

// @endcond


double zmin, zmax, delta;


ulong binStructure(const AtomicGroup& structure, Slices& slices) {
  ulong oob = 0;

  for (AtomicGroup::const_iterator a = structure.begin(); a != structure.end(); ++a) {
    double z = (*a)->coords().z();
    int bin = (z - zmin) * delta;
    if (bin < 0 || bin >= slices.size())
      ++oob;
    else {
      pAtom b(*a);
      b->coords().z(0.0);
      slices[bin].append(b);
    }
  }
  return(oob);
}



int main(int argc, char *argv[]) {

  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicSelection* sopts = new opts::BasicSelection("backbone");
  opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
  ToolOptions* topts = new ToolOptions;
  
  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(tropts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  AtomicGroup model = tropts->model;
  AtomicGroup subset = selectAtoms(tropts->model, sopts->selection);

  zmax = topts->zmax;
  zmin = topts->zmin;
  uint nbins = topts->nbins;
  delta = nbins / (zmax - zmin);
 
  vector<uint> frames = tropts->frameList();
  uint n = frames.size();

  
  RealMatrix M(nbins, n);
  
  uint col = 0;
  ulong out_of_bounds = 0;

  vector<double> avgs(nbins, 0.0);
  
  for (uint i=0; i<frames.size(); ++i) {
    tropts->trajectory->readFrame(frames[i]);
    tropts->trajectory->updateGroupCoords(subset);

    Slices slices(nbins);
    out_of_bounds += binStructure(subset, slices);
    for (uint j=0; j<nbins; ++j) {
      double d = 0.0;
      if (! slices[j].empty())
        d = topts->rgyr_mode ? slices[j].radiusOfGyration() : slices[j].radius();
      avgs[j] += d;
      M(j, col) = d;
    }
    ++col;
  }

  out_of_bounds /= n;
  
  for (uint i=0; i<nbins; ++i)
    avgs[i] /= n;

  vector<double> devs(nbins, 0.0);
  for (uint i=0; i<n; ++i)
    for (uint j=0; j<nbins; ++j) {
      double d = M(j, i) - avgs[j];
      devs[j] += d*d;
    }

  for (uint j=0; j<nbins; ++j)
    devs[j] = sqrt(devs[j]/(n-1));

  if (! topts->timeseries_filename.empty())
    writeAsciiMatrix(topts->timeseries_filename, M, hdr);

  cout << "# " << hdr << endl;
  cout << "# out of bounds = " << out_of_bounds << endl;
  cout << "# bin\tz\tavg\tstd\n";
  for (uint j=0; j<nbins; ++j) {
    double z = (j+0.5) * (zmax - zmin) / nbins + zmin;
    cout << boost::format("%d\t%f\t%f\t%f\n")
      % j
      % z
      % avgs[j]
      % devs[j];
  }
}
