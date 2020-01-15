/*
  dihedrals

  Computes the dihedral angle between each set of four atoms specified.
*/

/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008-2020 Tod D. Romo & Alan Grossfield
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

#include <boost/format.hpp>
#include <fstream>
#include <iostream>
#include <loos.hpp>
#include <regex>
#include <string>

using namespace std;
using namespace loos;
namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

typedef vector<AtomicGroup> vGroup;

// @cond TOOL_INTERNAL

string fullHelpMessage() {
  string s = "XXX";
  return s;
}

// these determine where the string containing the dihedral selections is split
const string quartet_delim = "|";
const string atom_delim = ",";

// C++ 11 regex split
// https://stackoverflow.com/questions/9435385/split-a-string-using-c11
vector<string> split(const string &input, const string &regex) {
  // passing -1 as the submatch index parameter performs splitting
  regex re(regex);
  sregex_token_iterator first{input.begin(), input.end(), re, -1}, last;
  return {first, last};
}

// split of split strings
vector<vector<string>> deep_split(const string &input,
                                  const string &outer_regex,
                                  const string &inner_regex) {
  vector<vector<string>> vector_of_frags;
  vector<string> outer_frags = split(input, outer_regex);
  // split the outer fragments and push these vectors into return object
  for (string &frag : outer_frags) {
    vector_of_frags.push_back(move(split(frag, inner_regex)));
  }
  return vector_of_frags;
}
class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions()
      : dihedral_sel_strings(""), pdb(""), tags(""), dihedral_sels{} {};
  // clang-format off
  void addGeneric(po::options_description &o) {
    o.add_options()
    ("dihedral-sel-strings,D", po::value<string>(&dihedral_sel_strings)->default_value(""),
     "Ordered quartets of selection strings; each quartet is delimited by '" 
     + quartet_delim + "', and each string within by '" 
     + atom_delim + "'.")
    ("pdb", po::value<string>(&pdb)->default_value(""),
     "Prefix to write PDBs for each dihedral selected from frame 1 of provided multi-traj.")
    ("tags,T", po::value<string>(&tags)->default_value(""), 
     "String of tags for each class of dihedral, separated by a '" + atom_delim + "'.");
  }
  // clang-format on

  string print() const {
    ostringstream oss;
    oss << boost::format("dihedral-sel-strings=%s,pdb=%s,tags=%s") %
               dihedral_sel_strings % pdb % tags;
  }

  bool postConditions(po::variables_map &map) {
    dihedral_sels = deep_split(dihedral_sel_strings, quartet_delim, atom_delim);
    for (auto d : dihedral_sels) {
      if (d.size() != 4) {
        throw(LOOSError("The following selection did not split to a quartet of "
                        "selections:\n"));
        for (auto s : d)
          cout cout << s << "\t";
        cout << "\n";
        return false;
      }
    }
    return true;
  }

  vector<vector<string>> dihedral_sels;
  string dihedral_sel_strings;
  string pdb;
  string tags;
};

// takes an atomic group for scope, and a vector of vectors of sel-strings.
// Corrects order of discovery  of each dihedral, and returns atomic group of
// dihedrals.
vector<vGroup> sels_to_dihedralAGs(const vector<vector<string>> &dihedral_sels,
                                   const AtomicGroup &scope) {
  // append to this for return later.
  vGroup dihedralAGs;
  for (auto dSels : dihedral_sels) {
    // first get a set of AGs that have all the atoms of the dihedral in them
    // They are likely to be in the order of the selection matched first
    // i.e. all the matches for selection 1, then all for 2, and so forth.
    AtomicGroup outoforder_dihedralType;
    for (auto sel : dSels) {
      outoforder_dihedralType += selectAtoms(scope, sel);
    }
    // This separates all non-connected atoms into separate atomic groups.
    vGroup dihedralTypeVector = outoforder_dihedralType.splitByMolecule();
    // reorder them here to match that provided by user
    // it may turn out this is unnecessary,
    // but the return order of selectAtoms calls is not specified.
    for (auto oo_D : dihedralTypeVector) {
      AtomicGroup reordered;
      for (auto sel : dSels) {
        reordered += SelectAtoms(oo_D, sel);
      }
      oo_D = move(reordered);
    }
    dihedralAGs.push_back(move(dihedralTypeVector));
  }
  return dihedralAGs;
}

int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);

  opts::BasicOptions *bopts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicSelection *sopts = new opts::BasicSelection("all");
  opts::MultiTrajOptions *mtopts = new opts::MultiTrajOptions;
  ToolOptions *topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(mtopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  // set up system for looping. Load coords from frame 0 into scope.
  AtomicGroup model = mtopts->model;
  AtomicGroup scope = selectAtoms(model, sopts->selection);
  pTraj traj = mtopts->trajectory;
  traj->updateGroupCoords(model);

  // figure out what dihedrals to track
  vector<vGroup> dihedrals = sels_to_dihedralAGs(topts->dihedral_sels, scope);

  // if tags provided, split those into vector.
  if (topts->tags) 
    vector<string> vtags = split(topts->tags, atom_delim)
  

  // if verbosity, and no pdbs were requested, then print each atomic group
  // found for each atom in each dihedral to stderr.
  if (bopts->verbosity > 0) {
    if (topts->pdb.empty()) {
      cerr << header;
      cerr << "# Following are the tab-delimited dihedral class selection "
              "strings and the atomic groups each produced:\n";
      for (uint i = 0; i < topts->dihedral_sels.size(); i++) {
        cerr << i << "\t";
        for (auto sel : topts->dihedral_sels.at(i)) {
          cerr << sel << "\t";
        }
        cerr << "\n[";
        for (auto ag : dihedrals.at(i)) {
          cerr << ag << ",";
        }
        cerr << "\b]\n";
      }
    }
  }

  // if PDB name string was given, write PDBs to indexed files by that prefix
  if (!topts->pdb.empty()) {
    for (uint i = 0; i < dihedrals.size(); i++) {
      for (uint j = 0; j < dihedrals.at(i).size(); i++) {
        PDB pdb = PDB::fromAtomicGroup(dihedrals[i, j]);
        pdb.remarks().add(to_string(j) + " from: " + dihedral_sels.at(i));
        ofstream pdbFile;
        pdbFile.open(topts->pdb + "_" + to_string(i) + "_" + to_string(j) +
                     ".pdb");
        pdbFile << pdb;
        pdbFile.close();
      }
    }
    PDB scopePDB = PDB::fromAtomicGroup(scope);
    ofstream scopeFile;
    scopeFile.open(topts->pdb + "_scope.pdb");
    scope.remarks().add(header);
    scopeFile << scopePDB;
  }

  cout << header;
  cout << "#\t";
  for (auto selset : dihedral_sels) {
    for (auto sel : selset) {
      cout << sel << ",";
    }
  }
  // Trajectory Loop here.
  while (traj->readFrame()) {
    traj->updateGroupCoords(model);
    for (auto)
  }
}