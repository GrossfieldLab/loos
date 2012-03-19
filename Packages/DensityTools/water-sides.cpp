/*
  water-sides.cpp

  Desc:
    Classifies a water as being on one side of the membrane or the
    other or inside the membrane (1 = upper, 0 = inside, -1 = lower).
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
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
#include <limits>
#include <boost/format.hpp>
#include <DensityTools.hpp>


using namespace std;
using namespace loos;



// @cond TOOLS_INTERNAL

typedef std::pair<double,double> Range;
typedef Math::Matrix<int, Math::ColMajor> Matrix;

Range membrane(0.0, 0.0);
string model_name, traj_name, selection_string;

enum Location { UPPER = 1, MEMBRANE = 0, LOWER = -1 };


string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\n"
    "\tClassify waters as above, below, or inside a membrane (based on Z-coordinate)\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\twater-sides constructs a TxN matrix where T is the size of the trajectory (# of frames)\n"
    "and N is the number of water molecules.  Each element has a value of 1 (above membrane),\n"
    "0 (inside membrane), or -1 (below membrane).  The classification of the water is based\n"
    "solely on it's z-coordinate and the range specified on the command line.\n"
    "\n"
    "\nEXAMPLES\n"
    "\twater-sides foo.pdb foo.dcd -15 15\n"
    "This example uses the default water selection (\"name == 'OH2'\") and places the\n"
    "membrane at -15 <= Z <= 15\n"
    "\n"
    "\twater-sides --selection 'name == \"HOH\"' foo.pdb foo.dcd -25 20\n"
    "This example picks all atoms called \"HOH\" as waters and places the membrane\n"
    "at -25 <= Z <= 20\n"
    ;

  return(msg);
}


class WaterSidesOptions : public opts::OptionsPackage {
public:

  WaterSidesOptions() : lower_bounds(0.0), upper_bounds(0.0) { }

  void addHidden(po::options_description& options) {
    options.add_options()
      ("lower", po::value<double>(&lower_bounds), "Lower leaflet bounds")
      ("upper", po::value<double>(&upper_bounds), "Upper leaflet bounds");
  }

  void addPositional(po::positional_options_description& options) {
    options.add("lower", 1);
    options.add("upper", 1);
  }

  bool check(po::variables_map& map) {
    return(! (map.count("lower") && map.count("upper")) );
  }

  string help() const { return("membrane-lower-bounds membrane-upper-bounds"); }
  string print() const {
    ostringstream oss;
    oss << boost::format("lower=%f, upper=%f") % lower_bounds % upper_bounds;
    return(oss.str());
  }

  double lower_bounds, upper_bounds;
};

// @endcond


Range parseRange(const string& s) {
  double a, b;
  int i = sscanf(s.c_str(), "%lf:%lf", &a, &b);
  if (i != 2) {
    cerr << "Parse error with " << s << endl;
    exit(-1);
  }
  return(Range(a,b));
}




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  opts::BasicOptions *basic_opts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicSelection *basic_selection = new opts::BasicSelection("name == 'OH2'");
  opts::TrajectoryWithFrameIndices *basic_traj = new opts::TrajectoryWithFrameIndices;
  WaterSidesOptions *my_opts = new WaterSidesOptions;

  opts::AggregateOptions options;
  options.add(basic_opts).add(basic_selection).add(basic_traj).add(my_opts);
  if (!options.parse(argc, argv)) {
    exit(0);
  }

  AtomicGroup model = basic_traj->model;
  pTraj traj = basic_traj->trajectory;
  AtomicGroup subset = selectAtoms(model, basic_selection->selection);
  vector<uint> frames = basic_traj->frameList();

  uint n = subset.size();
  uint m = frames.size();
  
  Matrix M(m,n+1);
  for (uint j=0; j<m; ++j) {
    traj->readFrame(frames[j]);
    traj->updateGroupCoords(subset);
    M(j, 0) = frames[j];
    for (uint i=0; i<n; ++i) {
      GCoord c = subset[i]->coords();
      Location l;
      if (c[2] > membrane.second)
        l = UPPER;
      else if (c[2] >= membrane.first)
        l = MEMBRANE;
      else
        l = LOWER;
      M(j, i+1) = l;
    }
  }

  writeAsciiMatrix(cout, M, hdr);
}
