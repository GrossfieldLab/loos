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
#include <boost/algorithm/string/replace.hpp>
#include <fstream>
#include <iostream>
#include <loos.hpp>
#include <regex>
#include <string>

using namespace std;
using namespace loos;
namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

// @cond TOOL_INTERNAL

string fullHelpMessage() {
  string s = "XXX";
  return s;
}

// these determine where the string containing the dihedral selections is split
const string quartet_delim = ":";
const string atom_delim = ",";
const string tag_delim = "_";
const string fsuffix = ".out";

// C++ 11 regex split
// https://stackoverflow.com/questions/9435385/split-a-string-using-c11
vector<string> split(const string &input, const string &regular_expression) {
  // passing -1 as the submatch index parameter performs splitting
  regex re(regular_expression);
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
    vector_of_frags.push_back(split(frag, inner_regex));
  }
  return vector_of_frags;
}
class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions()
      : dihedral_sels{}, dihedral_sel_strings(""), pdb(""), tags(""),
        prefix("dihedral"), quotes("p") {};
  // clang-format off
  void addGeneric(po::options_description& o) {
    o.add_options()
    ("dihedral-sel-strings,D", po::value<string>(&dihedral_sel_strings)->default_value(""),
     ("Ordered quartets of selection strings; each quartet is delimited by '"
     + quartet_delim + "', and each string within by '" 
     + atom_delim + "'.").c_str())
    ("pdb", po::value<string>(&pdb)->default_value(""),
     "Prefix to write PDBs for each dihedral selected from frame 1 of provided multi-traj.")
    ("tags,T", po::value<string>(&tags)->default_value(""), 
     ("String of tags for each class of dihedral, separated by a '" + atom_delim + "'.").c_str())
    ("prefix,p", po::value<string>(&prefix)->default_value("dihedral"),
     "Prefix for file names for each monitored dihedral."),
    ("swap-single-quotes,Q", po::value<string>(&quotes)->default_value("p"),
     "Swap single quote character in tags for some alternative. Provide single quote if no change desired..")
    ;
  }
  // clang-format on

  string print() const {
    ostringstream oss;
    oss << boost::format("dihedral-sel-strings=%s,pdb=%s,tags=%s,prefix=%s,quotes=%s") %
               dihedral_sel_strings % pdb % tags % prefix % quotes;
    return (oss.str());
  }

  bool postConditions(po::variables_map &map) {
    dihedral_sels = deep_split(dihedral_sel_strings, quartet_delim, atom_delim);
    for (auto d : dihedral_sels) {
      if (d.size() != 4) {
        throw(LOOSError("The following selection did not split to a quartet of "
                        "selections:\n"));
        for (auto s : d)
          cout << s << "\t";
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
  string prefix;
  string quotes;
};

// takes an atomic group for scope, and a vector of vectors of sel-strings.
// Corrects order of discovery  of each dihedral, and returns atomic group of
// dihedrals.
vector<vector<AtomicGroup>>
sels_to_dihedralAGs(const vector<vector<string>> &dihedral_sels,
                    const AtomicGroup &scope, const int verbosity) {

  // pick whether to puke up atomic group when group has wrong num elts.
  bool (*chkDihedralSize)(AtomicGroup &, vector<string>);
  if (verbosity > 0) {
    chkDihedralSize = [](AtomicGroup &oo_D, vector<string> sels) -> bool {
      if (oo_D.size() != 4) {
        cerr << "WARNING: dihedral specification found " << oo_D.size();
        cerr << " atoms, not 4 in selection string set: \n\t";
        for (auto sel : sels)
          cerr << sel << ", ";
        cerr << "\b\b\n";
        cerr << "Offending group: \n";
        cerr << oo_D;
        cerr << "\nDROPPING THIS GROUP AND PROCEEDING.\n";
        return true;
      } else {
        AtomicGroup reordered;
        for (auto sel : sels)
          reordered += selectAtoms(oo_D, sel);

        oo_D = move(reordered);
        cerr << "included group of size: " << to_string(reordered.size())
             << "\n";
        return false;
      }
    };
  } else {
    chkDihedralSize = [](AtomicGroup &oo_D, vector<string> sels) -> bool {
      if (oo_D.size() != 4)
        return true;
      else {
        AtomicGroup reordered;
        for (auto sel : sels)
          reordered += selectAtoms(oo_D, sel);

        oo_D = move(reordered);
        return false;
      }
    };
  }

  // append to this for return later.
  vector<vector<AtomicGroup>> dihedralAGs;
  for (auto dSels : dihedral_sels) {
    // first get a set of AGs that have all the atoms of the dihedral in them
    // They are likely to be in the order of the selection matched first
    // i.e. all the matches for selection 1, then all for 2, and so forth.
    AtomicGroup outoforder_dihedralType;
    for (auto sel : dSels) {
      outoforder_dihedralType += selectAtoms(scope, sel);
    }
    // This separates all non-connected atoms into separate atomic groups.
    vector<AtomicGroup> dihedralInstances =
        outoforder_dihedralType.splitByMolecule();
    // reorder them here to match that provided by user
    // it may turn out this is unnecessary,
    // but the return order of selectAtoms calls is not specified.
    // Remove any AGs that didn't manage to contain four atoms after the
    // split.
    dihedralInstances.erase(remove_if(dihedralInstances.begin(),
                                      dihedralInstances.end(),
                                      [&](AtomicGroup &oo_D) -> bool {
                                        return (*chkDihedralSize)(oo_D, dSels);
                                      }),
                            dihedralInstances.end());
    dihedralAGs.push_back(move(dihedralInstances));
  }
  return dihedralAGs;
}

int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);

  opts::BasicOptions *bopts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicSelection *sopts =
      new opts::BasicSelection("backbone && !hydrogen");
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
  vector<vector<AtomicGroup>> dihedrals =
      sels_to_dihedralAGs(topts->dihedral_sels, scope, bopts->verbosity);

  // make file names, either from scratch or by adding to user appended tags.
  vector<vector<shared_ptr<ofstream>>> vv_filePtrs;
 // if user supplied tags for file names, use those with reduced dihedral name info.
  if (topts->tags.empty()) {
    int resid;
    for (auto dihedralType : dihedrals) {
      vector<shared_ptr<ofstream>> v_filePtrs;
      for (auto dihedral : dihedralType) {
        resid = dihedral[0]->resid();
        string tag;
        for (auto patom : dihedral) {
          // put a residue number with the name for each atom not from residue
          // of atom zero.
          string name = patom->name();
          boost::replace_all(name, "\'", topts->quotes);
          if (resid != patom->resid())
            tag = tag_delim + to_string(patom->resid()) + name;
          else
            tag = tag_delim + name;
        }
        auto p_ofstream =
            make_shared<ofstream>(topts->prefix + tag_delim + tag + fsuffix);
        *(p_ofstream) << "# " << header << "\n";
        v_filePtrs.push_back(p_ofstream);
      }
      vv_filePtrs.push_back(move(v_filePtrs));
    }
  } else {
    vector<string> user_tags = split(topts->tags, atom_delim);
    for (uint i = 0; i < user_tags.size(); i++) {
      vector<shared_ptr<ofstream>> v_filePtrs;
      for (auto dihedral : dihedrals.at(i)) {
        string tag = user_tags.at(i);
        tag += tag_delim + to_string(dihedral[0]->resid());
        for (auto patom : dihedral){
          string name = patom->name();
          boost::replace_all(name, "\'", topts->quotes);
          tag += tag_delim + name; // append atom names to tag with
        }
                                            // tag delimiter
        auto p_ofstream =
            make_shared<ofstream>(topts->prefix + tag_delim + tag + fsuffix);
        *p_ofstream << "# " << header << "\n";
        v_filePtrs.push_back(p_ofstream);
      }
      vv_filePtrs.push_back(move(v_filePtrs));
    }
  }

  // if verbosity, and no pdbs were requested, then print each atomic group
  // found for each atom in each dihedral to stderr.
  if (bopts->verbosity > 0) {
    if (topts->pdb.empty()) {
      cerr << "# " << header << "\n";
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
      for (uint j = 0; j < (dihedrals.at(i)).size(); j++) {
        PDB pdb = PDB::fromAtomicGroup(dihedrals.at(i).at(j));
        string rmks = to_string(j) + " from: ";
        for (auto sel : topts->dihedral_sels.at(i))
          rmks += sel + ", ";

        rmks += "\b\b";

        pdb.remarks().add(rmks);

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
    scopePDB.remarks().add(header);
    scopeFile << scopePDB;
  }

  // Trajectory Loop here.
  double dihedral_angle;
  while (traj->readFrame()) {
    traj->updateGroupCoords(model);
    for (uint typeIdx = 0; typeIdx < dihedrals.size(); typeIdx++) {
      for (uint dIdx = 0; dIdx < dihedrals.at(typeIdx).size(); dIdx++) {
        dihedral_angle = Math::torsion(dihedrals[typeIdx][dIdx][0]->coords(),
                                       dihedrals[typeIdx][dIdx][1]->coords(),
                                       dihedrals[typeIdx][dIdx][2]->coords(),
                                       dihedrals[typeIdx][dIdx][3]->coords());
        *(vv_filePtrs[typeIdx][dIdx])
            << traj->currentFrame() << "\t" << dihedral_angle << "\n";
      }
    }
  }
  // close all these output files now that we're done looping over traj
  for (auto v_fps : vv_filePtrs)
    for (auto p_ofs : v_fps)
      p_ofs->close();
}