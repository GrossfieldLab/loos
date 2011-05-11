/*
  average.cpp

  Computes the average structure post-aligning...
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
#include <boost/program_options.hpp>

using namespace std;
namespace opts = loos::OptionsFramework;
using namespace loos;





// @cond TOOL_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions(const string& s) : avg_string(s) { }

  void addGeneric(opts::po::options_description& o) {
    o.add_options()
      ("average", opts::po::value<string>(&avg_string)->default_value(avg_string), "Average over this selection");
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

  opts::BasicOptions* bopts = new opts::BasicOptions;
  opts::BasicSelectionOptions* sopts = new opts::BasicSelectionOptions("");
  opts::BasicTrajectoryOptions* trajopts = new opts::BasicTrajectoryOptions;
  ToolOptions* toolopts = new ToolOptions("!hydrogen || segid == 'SOLV' || segid == 'BULK'");

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(trajopts).add(toolopts);
  if (!options.parse(argc, argv))
    exit(-1);



  AtomicGroup model = createSystem(trajopts->model_name);

  AtomicGroup avg_subset = selectAtoms(model, toolopts->avg_string);
  cerr << "Averaging over " << avg_subset.size() << " atoms.\n";

  pTraj traj = createTrajectory(trajopts->traj_name, model);
  vector<uint> indices = opts::assignFrameIndices(traj, trajopts->frame_index_spec, trajopts->skip);
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
