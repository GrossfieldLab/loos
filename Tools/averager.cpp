/*
  averager.cpp

  Computes the average structure post-aligning

  Usage:
    averager [options] model traj >average.pdb


  The --selection option determines which atoms are used for the
  iterative alignment.  The --average option determines which atoms
  are actually averaged and written out.
*/



/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo
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
    "Compute an average structure from a trajectory\n"
    "\n"
    "DESCRIPTION\n"
    "\taverager writes out a PDB for the average structure from a trajectory.  If a selection\n"
    "is given (--selection), then the trajectory is first iteratively aligned to an optimal\n"
    "average structure (see aligner).  The '--average' option takes an optional selection that\n"
    "defines what atoms are averaged and written out, otherwise all non-hydrogen and non-solvent\n"
    "atoms are used.  Note that solvent is selected by a segid of either 'BULK' or 'SOLVENT'.\n"
    "If your system uses a different identifier, you will want to explicitly give a selection\n"
    "for the --average option\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\taverager model.pdb traj.dcd >average.pdb\n"
    "This assumes the trajectory is already aligned and puts the average structure in average.pdb\n"
    "Hydrogens and solvent atoms are excluded.\n"
    "\n"
    "\taverager --selection 'name == \"CA\"' model.pdb traj.dcd >average.pdb\n"
    "Aligns the trajectory first using all alpha-carbons\n"
    "\n"
    "\taverager --selection 'name == \"CA\"' --average 'resid <= 20' model.pdb traj.dcd >average.pdb\n"
    "Aligns the trajectory using alpha-carbons, but only averages the first 20 residues and outputs\n"
    "them to average.pdb\n"
    "\n"
    "SEE ALSO\n"
    "\taligner\n";

  return(msg);
}


class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions(const string& s) : avg_string(s) { }

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("average", po::value<string>(&avg_string)->default_value(avg_string), "Average over this selection");
  }

  string print() const {
    ostringstream oss;

    oss << boost::format("avg_string='%s'") % avg_string;
    return(oss.str());
  }


  string avg_string;
};

// @endcond





int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicSelection* sopts = new opts::BasicSelection("");
  opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
  ToolOptions* toolopts = new ToolOptions("!(hydrogen || segid == 'SOLV' || segid == 'BULK')");

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(tropts).add(toolopts);
  if (!options.parse(argc, argv))
    exit(-1);



  AtomicGroup model = tropts->model;

  AtomicGroup avg_subset = selectAtoms(model, toolopts->avg_string);
  cerr << "Averaging over " << avg_subset.size() << " atoms.\n";

  pTraj traj = tropts->trajectory;
  vector<uint> indices = tropts->frameList();
  cerr << "Using " << indices.size() << " frames from the trajectory...\n";

  // First, align...
  vector<XForm> xforms;
  if (sopts->selection.empty()) {

    cerr << "Skipping alignment...\n";
    for (uint i=0; i<indices.size(); ++i)
      xforms.push_back(XForm());

  } else {

    AtomicGroup align_subset = selectAtoms(model, sopts->selection);
    cerr << "Aligning with " << align_subset.size() << " atoms.\n";

    boost::tuple<vector<XForm>, greal, int> result = iterativeAlignment(align_subset, traj, indices);
    xforms = boost::get<0>(result);
    double rmsd = boost::get<1>(result);
    int niters = boost::get<2>(result);
    cerr << boost::format("Aligned in %d iterations with final error of %g.\n") % niters % rmsd;
  }

  // Now re-read the average subset
  cerr << "Averaging...\n";
  AtomicGroup avg = averageStructure(avg_subset, xforms, traj, indices);
  PDB avgpdb = PDB::fromAtomicGroup(avg);
  avgpdb.pruneBonds();
  avgpdb.remarks().add(header);
  cout << avgpdb;
}
