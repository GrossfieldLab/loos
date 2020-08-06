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

#include <boost/algorithm/string/replace.hpp>
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

// @cond TOOL_INTERNAL

// these determine where the string containing the dihedral selections is split
const string quartet_delim = ":";
const string atom_delim = ",";
const string tag_delim = "_";
const string fsuffix = ".out";

// clang-format off
const string msg = 
"SYNOPSIS\n"
"\n"
"This tool is designed to allow the tracking of classes of dihedral angles \n"
"specified by atom selection. \n"
"\nDESCRIPTION\n"
"\n"
"Unlike the torsion tool, also in LOOS, this tool is designed to track the\n"
"dihedral angle between chemically connected groups of four atoms. The \n"
"original intention for the tool was to monitor classes of customarily\n"
"defined dihedrals that might exist in a large number of residues in one pass.\n"
"For example, one could use this tool to monitor all of the phi and psi\n"
"backbone dihedrals in a protein, making only one pass through the trajectory as\n"
"one did so. The tool creates a file name for each dihedral angle chosen for \n"
"monitoring, and writes the frame number and the angle out in two columns, \n"
"separated by white space, for each frame provided to the tool. How these names \n"
"are created, how many classes of dihedral to monitor, and what frames to \n"
"consider from the input trajectory(ies) are all configurable. Because it \n"
"handles output through a number of out files, it can be wise to create \n"
"a subdirectory that will contain the mess.\n"
" \n"
"The --selection flag controls the scope of the search for dihedrals to monitor.\n"
" So in the aforementioned protein example, if you only wanted to monitor the \n"
"phi of the first five residues of some protein, you would provide a selection \n"
"string like 'resid < 6' (assuming of course that your protein's residues are \n"
"the first such in the overall list of residues, which is commonly the case).  \n"
" \n"
"Several of the flags are from LOOS classes devoted to providing basic tool \n"
"functionality, and they work the same as in other tools. For example, \n"
"trajectories are read using a MultiTrajectory, and so the skip, stride, and \n"
"range flags all do what they do for multi-trajectory based tools. This is also \n"
"why you can provide an arbitrary number of trajectories to this tool, and it \n"
"will gracefully treat them as one long trajectory. \n"
" \n"
"The --dihedral-sel-strings flag is obligate. It should be a string that \n"
"provides a list of atom selections in quartets separated by a '"+atom_delim+"'. Each \n"
"selection string should grab only one atom so that each quartet selects four \n"
"atoms, in the order that you would like them fed to the loos::Math::torsion() \n"
"function. If you'd like to monitor multiple types of dihedral, even if it's the\n"
" same dihedral across different residues (for example, chi, the glycosidic \n"
"dihedral in nucleic acids) you can include multiple quartets by interspersing \n"
"'"+quartet_delim+"' between each quartet. For example, to select the chi dihedral in nucleic \n"
"acids you could write: \n"
" \n"
"    --dihedral-sel-strings $\'name == \"O4'\""+atom_delim+"  name == \"C1'\""+atom_delim+"  name == \"N9\""+atom_delim+"  \\\n"
"name == \"C4\" "+quartet_delim+" name == \"O4'\""+atom_delim+"  name == \"C1'\""+atom_delim+"  name == \"N1\""+atom_delim+"  name == \"C2\"\' \n"
" \n"
"Noting that the four selection strings before the '"+quartet_delim+"' are for purine chis, and \n"
"the four after are for pyrimidine chis. In the case of nucleic acids, which \n"
"usually have the \"\'\" character in the atom name, it can be very helpful to \n"
"put the arguments to this tool in a config file. See the LOOS online docs for \n"
"how to go about that. \n"
" \n"
"The --pdb flag is for debugging. If you want to use it, provide a prefix by \n"
"which to name the reported pdb files. It takes the first frame of the multi-\n"
"trajectory and writes out the scope, and each four atom sequence it found as \n"
"separate PDB files, prefixed with the provided argument. For each PDB created \n"
"thus, it numbers the files first by dihedral class, then by which element in \n"
"the class it is. So if you provide the 'test' as an argument, your PDBs might \n"
"look like: \n"
" \n"
"    test_x_y.pdb \n"
" \n"
"Where the contents will be the yth dihedral of type x found. To get a nice \n"
"visual representation of how the selection went, I like to say 'pymol *.pdb' in\n"
" the subdirectory I made for this analysis, then show all as 'sticks/licorice',\n"
" and overwrite that setting for just the scope with 'lines'. This makes it \n"
"patently clear where the dihedrals being tracked will be in the molecule. \n"
" \n"
"The --tags option is for providing tags that correspond to each class of \n"
"dihedrals monitored by each quartet. Each tag provides an infix name that \n"
"corresponds to the selection string that is in that position in the '"+atom_delim+"' \n"
"separated list of --dihedral-sel-strings. For the chi example: \n"
" \n"
"    --tags 'chi_R,chi_Y' \n"
" \n"
"Since the first of the two quartets corresponds to purines and the second to \n"
"pyrimidines. If you do provide this argument, it needs to have the same number \n"
"of ',' separated strings as you've provided quartets above. If you elect not to\n"
" provide it, then a tag is fabricated from the residue name, resid, and each of\n"
" the atoms selected, separated by '"+tag_delim+"'. If some of the atoms cross into another \n"
"residue, those atom names will have the resid of that neighboring residue \n"
"appearing after the name. If you do provide tags, then the tag, followed by the\n"
" resid of the first atom, then the names of the atoms in that particular \n"
"dihedral will be the filename instead, also separated by '"+tag_delim+"'. In the chi \n"
"example, because of the tags provided, an output file might look like the \n"
"following: \n"
" \n"
"    roc_chi_R_1_O4p_C1p_N9_C4" +fsuffix+" \n"
" \n"
"Note that the primes have been replaced by the letter p, which can be changed \n"
"(even back to _shudder_ a \') if the user specifies the --swap-single-quotes \n"
"flag. \n"
" \n"
"The --prefix flag is a string that precedes all the dihedral time series file \n"
"names (aside from the output caused by --pdbs) This permits exclusive names for\n"
" different runs of the program and helps keep things organized. I often use a \n"
"system specifying prefix.\n"
"\n"
"EXAMPLE\n"
"\n"
"dihedrals \\\n--dihedral-sel-strings $\'name == \"O4'\""+atom_delim+"  name == \"C1'\""+atom_delim+"  name == \"N9\""+atom_delim+" \\\n"
"name == \"C4\" "+quartet_delim+" name == \"O4'\""+atom_delim+"  name == \"C1'\""+atom_delim+"  name == \"N1\""+atom_delim+"  name == \"C2\"\' \\\n"
"--tags  'chi_Y,chi_R' --selection 'resid < 6' --prefix nucX nuc.pdb nuc.dcd\n"
"\n"
"This should do the calculation discussed in the description above. In particular\n"
"it will look for dihedrals matching the conventional names for chi from \n"
"purines and pyrimidines, writing each instance of these classes out to different\n"
"output files with names based on --prefix.\n"
"\n"
"POTENTIAL COMPLICATIONS\n"
"\n"
"Verbosity and the --pdb flag help diagnose problems with dihedral selections.\n"
"This is a very good thing to check with all tools, but especially here, where \n"
"results could look right but be wrong with selection strings that are subtly off.\n"
"\n"
"Another thing to bear in mind is that the model needs connectivity. One can\n"
"remedy this with the --infer-connectivity flag, but use caution. That inference\n"
"can be low quality if one gets unlucky with the first frame in the file, since\n"
"it is based on how far apart atoms are from one another. Regardless of what\n"
"is provided for this flag, if connectivity information is found then none will\n"
"be inferrd.\n"
;
// clang-format on

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
        prefix("dihedral"), quotes("p"){};
  // clang-format off
  void addGeneric(po::options_description& o) {
    o.add_options()
    ("dihedral-sel-strings,D", po::value<string>(&dihedral_sel_strings)->default_value(""),
     ("Ordered quartets of selection strings; each quartet is delimited by '"
     + quartet_delim + "', and each string within by '" 
     + atom_delim + "'.").c_str())
    ("infer-connectivity", po::value<float>(&bondlength)->default_value(-1), 
    "Infer connectivity using provided distance for models lacking this. ALERT: uses hard distance cutoff on first frame of traj to infer connectivity. Only does this for values greater than zero.")
    ("pdb", po::value<string>(&pdb)->default_value(""),
     "Prefix to write PDBs for each dihedral selected from frame 1 of provided multi-traj.")
    ("tags,T", po::value<string>(&tags)->default_value(""), 
     ("String of tags for each class of dihedral, separated by a '" + atom_delim + "'.").c_str())
    ("prefix,p", po::value<string>(&prefix)->default_value("dihedral"),
     "Prefix for file names for each monitored dihedral.")
    
    ("swap-single-quotes,Q", po::value<string>(&quotes)->default_value("p"),
     "Swap single quote character in outfile names for some alternative. Provide single quote if no change desired.");
  }
  // clang-format on

  string print() const {
    ostringstream oss;
    oss << boost::format("dihedral-sel-strings=%s,pdb=%s,tags=%s,prefix=%s,"
                         "quotes=%s,bondlength=%d") %
               dihedral_sel_strings % pdb % tags % prefix % quotes % bondlength;
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
  float bondlength;
};

// takes an atomic group for scope, and a vector of vectors of sel-strings.
// Corrects order of discovery  of each dihedral, and returns atomic group of
// dihedrals.
vector<vector<AtomicGroup>>
sels_to_dihedralAGs(const vector<vector<string>> &dihedral_sels,
                    const AtomicGroup &scope, const int verbosity) {

  // pick whether to puke up atomic group when group has wrong num elts.
  bool (*chkSizeReorder)(AtomicGroup &, vector<string> &);

  if (verbosity > 0) {
    chkSizeReorder = [](AtomicGroup &oo_D, vector<string> &sels) -> bool {
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
      } else { // append AG, correctly reordered.
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
    chkSizeReorder = [](AtomicGroup &oo_D, vector<string> &sels) -> bool {
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
    dihedralInstances.erase(
        remove_if(dihedralInstances.begin(), dihedralInstances.end(),
                  [&dSels, chkSizeReorder](AtomicGroup &oo_D) -> bool {
                    return (*chkSizeReorder)(oo_D, dSels);
                  }),
        dihedralInstances.end());
    dihedralAGs.push_back(move(dihedralInstances));
  }
  return dihedralAGs;
}

int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);

  opts::BasicOptions *bopts = new opts::BasicOptions(msg);
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
  if (model.hasBonds()) {
  } else if (topts->bondlength > 0)
    model.findBonds(topts->bondlength);
  else
    throw(LOOSError(
        "Model does not appear to have chemical connectivity, and "
        "infer-connectivity has not been set to a positive value.\n"));
  AtomicGroup scope = selectAtoms(model, sopts->selection);
  pTraj traj = mtopts->trajectory;
  traj->updateGroupCoords(model);

  // figure out what dihedrals to track
  vector<vector<AtomicGroup>> dihedrals =
      sels_to_dihedralAGs(topts->dihedral_sels, scope, bopts->verbosity);

  // make file names, either from scratch or by adding to user appended tags.
  vector<vector<shared_ptr<ofstream>>> vv_filePtrs;
  // if user supplied tags for file names, use those with reduced dihedral name
  // info.
  if (topts->tags.empty()) {
    int resid;
    for (auto dihedralType : dihedrals) {
      vector<shared_ptr<ofstream>> v_filePtrs;
      for (auto dihedral : dihedralType) {
        resid = dihedral[0]->resid();
        string tag(dihedral[0]->resname() + to_string(resid));
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
        for (auto patom : dihedral) {
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
