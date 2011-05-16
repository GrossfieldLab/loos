/*
  perturb-structure.cpp

  Apply a random perturbation to a structure (random directions, fixed magnitude)
*/



/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2010, Tod D. Romo
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
#include <boost/format.hpp>
#include <boost/program_options.hpp>


using namespace std;
using namespace loos;
namespace opts = loos::OptionsFramework;



uint seed;
string selection;
string model_name;
double magnitude;


// @cond TOOL_INTERNAL

class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() : seed(0), magnitude(0.0) { }

  void addGeneric(opts::po::options_description& o) {
    o.add_options()
      ("seed", opts::po::value<uint>(&seed)->default_value(seed), "Random number seed (0 = use current time)");
  }
  
  void addHidden(opts::po::options_description& o) {
    o.add_options()
      ("magnitude", opts::po::value<double>(&magnitude), "magnitude");
  }

  void addPositional(opts::po::positional_options_description& pos) {
    pos.add("magnitude", 1);
  }

  bool check(opts::po::variables_map& map) {
    return(!map.count("magnitude"));
  }

  string help() const { return("magnitude"); }
  string print() const {
    ostringstream oss;
    oss << boost::format("seed=%d, magnitude=%f") % seed % magnitude;
    return(oss.str());
  }

  bool postConditions(opts::po::variables_map& map) {
    if (seed == 0)
      randomSeedRNG();
    else {
      base_generator_type& rng = rng_singleton();
      rng.seed(seed);
    }

    return(true);
  }

  uint seed;
  double magnitude;
};

// @endcond






int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  opts::BasicOptions* bopts = new opts::BasicOptions;
  opts::BasicSelectionOptions* sopts = new opts::BasicSelectionOptions;
  opts::ModelWithCoordsOptions* mwcopts = new opts::ModelWithCoordsOptions;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(mwcopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  AtomicGroup model = opts::loadStructureWithCoords(mwcopts->model_name, mwcopts->coords_name);
  AtomicGroup subset = selectAtoms(model, sopts->selection);

  subset.perturbCoords(topts->magnitude);
  PDB pdb = PDB::fromAtomicGroup(model);
  pdb.remarks().add(hdr);

  cout << pdb;
}
