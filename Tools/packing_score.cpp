/*
  packing_score.cpp

  (c) 2013 Alan Grossfield 
           Department of Biochemistry and Biophysics
           University of Rochster School of Medicine and Dentistry


*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2011 Tod D. Romo
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



string sel1, sel2;
bool normalize = false;


// The following conditional prevents this class from appearing in the
// Doxygen documentation for LOOS:
//
// @cond TOOLS_INTERNAL

string fullHelpMessage(void) {
    string msg = 
"\n"
"SYNOPSIS\n"
"\n"
"Compute the packing score, a measure of contact between two selection.\n"
"\n"
"DESCRIPTION\n"
"\n"
"This tool computes the packing score, a simple measure of the contact between\n"
"two selections, over the course of a trajectory.  The packing score is defined\n"
"as the sum of the inverse of the pairwise distance raised to the sixth power,\n"
"computed over all pairs of atoms in 2 selections.  In essence, you can think \n"
"of it as the attractive component of the van der Waal's interaction with sigma \n"
"and epsilon set to 1.\n"
"The packing score was originally defined in Grossfield, A., Feller, S. E., \n"
"Pitman, M. C., A role for direct interactions in the modulation of rhodopsin \n"
"by omega-3 polyunsaturated lipids, Proc. Nat. Acad. Sci. USA, 2006, 103, \n"
"4888-4893\n"
"\n"
"Required Flags\n"
"   --sel1 'selection string' : the first selection of atoms\n"
"   --sel2 'selection string' : the second selection of atoms\n"
"\n"
"Options\n"
"\n"
"   --skip  N   : skip the first N frames from the trajectory\n"
"   --normalize : if specified, divides the packing score by the product of the \n"
"                 number of atoms in selection 1 and in selection 2\n"
"\n"
"EXAMPLE\n"
"\n"
"An example command line would be:\n"
"\n"
"packing_score --skip 20 --normalize --sel1 'resid >= 35 && resid <= 64 && segname == \"RHOD\"' --sel2 'resid>=71 && resid<=100 && segname == \"RHOD\"' rhod_namd.psf rhod_control.dcd\n"
"\n"
"which would compute the packing score between two chunks of the segment \n"
"\"RHOD\", one residues 35-64 and the other 71-100, normalizing the value,\n"
"and skipping the first 20 frames.\n"
"\n"
"HINTS\n"
"-- One trick to speed things up is to add \"&& !hydrogen\" to your \n"
"   selections.  The answers with and without the hydrogens should be\n"
"   almost perfectly proportional, but you'll greatly reduce the number\n"
"   of distance calculations\n"
"\n"
"-- The program verifies that the two selections don't have any atoms in \n"
"   common, and quits immediately if they do, since then the packing score\n"
"   would be infinite.\n"
"\n"
"-- The normalize option could be useful for identifying the pieces of an\n"
"   interface that pack tightly in a size-independent manner.  Otherwise,\n"
"   large chunks (e.g. tryptophans vs. alanines) will tend to have higher \n"
"   scores just because they have more atoms.\n"
"\n";
    return(msg);

}


class ToolOptions : public opts::OptionsPackage {
public:

  // Change these options to reflect what your tool needs
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("sel1", po::value<string>(&sel1), "selection 1")
      ("sel2", po::value<string>(&sel2), "selection 2")
      ("normalize", po::bool_switch()->default_value(false), "Normalize the score by the number of pairs")
      ;
  }

  // The print() function returns a string that describes what all the
  // options are set to (for logging purposes)
  string print() const {
    ostringstream oss;
    oss << boost::format("sel1=%s, sel2=%s") % sel1 % sel2;
    return(oss.str());
  }

  bool postConditions(po::variables_map &map) {
      if (sel1.empty() || sel2.empty()) {
          cerr << "Error: must specify --sel1 and --sel2" << endl;
          return(false);
      }
              
      normalize = map["normalize"].as<bool>();

    return(true);
  }

};
// @endcond
// ----------------------------------------------------------------



int main(int argc, char *argv[]) {
  
  string header = invocationHeader(argc, argv);
  
  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicTrajectory* tropts = new opts::BasicTrajectory;
  ToolOptions* topts = new ToolOptions;
  opts::AggregateOptions options;
  options.add(bopts).add(tropts).add(topts);

  // Parse the command-line.  If an error occurred, help will already
  // be displayed and it will return a FALSE value.
  if (!options.parse(argc, argv))
    exit(-1);

  cout << "# " << header << endl;

  AtomicGroup model = tropts->model;
  pTraj traj = tropts->trajectory;
  AtomicGroup set1 = selectAtoms(model, sel1);
  AtomicGroup set2 = selectAtoms(model, sel2);

  // Check that the two groups don't overlap, since you'll get
  // a NaN packing score.
  AtomicGroup overlap = set1.intersect(set2);
  if (overlap.size()) {
      cerr << "Error: the two selections have the following atoms in common:" 
           << endl;
      for (uint i=0; i < overlap.size(); i++) {
          cerr << *(overlap[i]) << endl;
      }
      exit(1);
  }

  uint num_pairs = 0;
  if (normalize) {
      num_pairs = set1.size() * set2.size();
  }
      

  // Now iterate over all frames in the trajectory (excluding the skip
  // region)
  uint frame = tropts->skip;
  while (traj->readFrame()) {
    traj->updateGroupCoords(tropts->model);

    double score = 0.0;

    for (uint i=0; i<set1.size(); i++) {
        for (uint j=0; j<set2.size(); j++) {
            GCoord diff = set1[i]->coords() - set2[j]->coords();
            double dist2 = diff.length2();
            double dist6 = dist2*dist2*dist2;
            score += 1.0/dist6;
        }
    }

    if (normalize) {
        score /= num_pairs;
    }
  
    cout << frame << "\t" << score << endl;
    frame++;
  }

}
