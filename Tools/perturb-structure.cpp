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

using namespace std;
using namespace loos;
namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;




// @cond TOOL_INTERNAL


string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "Randomly perturb atom coordinates in a model\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tThis tool randomly perturbs the coordinates in a model.  A subset may be selected\n"
    "and perturbed, in which case the entire model is still written out.\n"
    "\n"
    "NOTES\n"
    "\tRequires a model with coordinates\n"
    "\n";

  return(msg);
}



class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() : seed(0), magnitude(0.0) { }

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("seed", po::value<uint>(&seed)->default_value(seed), "Random number seed (0 = use current time)");
  }
  
  void addHidden(po::options_description& o) {
    o.add_options()
      ("magnitude", po::value<double>(&magnitude), "magnitude");
  }

  void addPositional(po::positional_options_description& pos) {
    pos.add("magnitude", 1);
  }

  bool check(po::variables_map& map) {
    return(!map.count("magnitude"));
  }

  string help() const { return("magnitude"); }
  string print() const {
    ostringstream oss;
    oss << boost::format("seed=%d, magnitude=%f") % seed % magnitude;
    return(oss.str());
  }

  bool postConditions(po::variables_map& map) {
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
  opts::BasicSelection* sopts = new opts::BasicSelection;
  opts::ModelWithCoords* mwcopts = new opts::ModelWithCoords;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(mwcopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  AtomicGroup model = mwcopts->model;
  AtomicGroup subset = selectAtoms(model, sopts->selection);

  subset.perturbCoords(topts->magnitude);
  PDB pdb = PDB::fromAtomicGroup(model);
  pdb.remarks().add(hdr);

  cout << pdb;
}
