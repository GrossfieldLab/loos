/*
  rmsfit.cpp

  Superimposes one structure upon another
*/



/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009 Tod D. Romo
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

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;


using namespace std;
using namespace loos;

// Globals...blech...

string source_selection, target_selection, apply_selection;


// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("apply,A", po::value<string>(&apply_selection)->default_value("all"), "Subset of source model to apply transformation to")
      ("source,S", po::value<string>(&source_selection)->default_value("name == 'CA'"), "Subset of the source model to align with")
      ("target,T", po::value<string>(&target_selection)->default_value("name == 'CA'"), "Subset of the target model to align with");

  }

  void addHidden(po::options_description& o) {
    o.add_options()
      ("source_name", po::value<string>(&source_name), "Source model filename")
      ("target_name", po::value<string>(&target_name), "Target model filename");
  }

  void addPositional(po::positional_options_description& p) {
    p.add("source_name", 1);
    p.add("target_name", 1);
  }

  bool check(po::variables_map& vm) {
    return(!(vm.count("source_name") && vm.count("target_name")));
  }

  bool postConditions(po::variables_map& vm) {
    source_model = createSystem(source_name);
    target_model = createSystem(target_name);

    return(true);
  }

  string help() const { return("source-filename target-filename"); }
  string print() const {
    ostringstream oss;
    oss << boost::format("apply='%s', source='%s', target='%s, source_name='%s', target_name='%s'")
      % apply_selection
      % source_selection
      % target_selection
      % source_name
      % target_name;
    return(oss.str());
  }


  string source_name;
  string target_name;

  AtomicGroup source_model;
  AtomicGroup target_model;
};

// @endcond



int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  AtomicGroup source_model = topts->source_model;
  AtomicGroup target_model = topts->target_model;

  AtomicGroup source_subset = selectAtoms(source_model, source_selection);
  AtomicGroup apply_subset = selectAtoms(source_model, apply_selection);
  if (apply_selection != "all")
    apply_subset.clearBonds();

  AtomicGroup target_subset = selectAtoms(target_model, target_selection);

  if (source_subset.size() != target_subset.size()) {
    cerr << boost::format("ERROR - The source subset has %d atoms but the target subset has %d atoms.  The MUST be equal")
      % source_subset.size()
      % target_subset.size();
    exit(-10);
  }

  GMatrix M = source_subset.superposition(target_subset);
  XForm W(M);
  apply_subset.applyTransform(W);
  if (bopts->verbosity) {
    double d = source_subset.rmsd(target_subset);
    cerr << "Final RMSD = " << d << endl;
  }
  
  PDB pdb = PDB::fromAtomicGroup(apply_subset);
  pdb.remarks().add(hdr);
  cout << pdb;
}
