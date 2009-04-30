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
namespace po = boost::program_options;
using namespace loos;



string align_string;
string avg_string;
double alignment_tol;
string model_name, traj_name;
vector<uint> indices;

void parseOptions(int argc, char *argv[]) {
  vector<string> ranges;

  try {
    po::options_description generic("Allowed options", 120);
    generic.add_options()
      ("help", "Produce this help message")
      ("align,a", po::value<string>(&align_string),"Align using this selection (or skip aligning)")
      ("average,A", po::value<string>(&avg_string)->default_value("!(hydrogen || segid == 'SOLV' || segid == 'BULK')"), "Average over this selection")
      ("range,r", po::value< vector<string> >(&ranges), "Range of frames to average over (Octave-style)");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value<string>(&traj_name), "Trajectory filename");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("traj", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("traj"))) {
      cerr << "Usage- averager [options] model-name trajectory-name\n";
      cerr << generic;
      exit(-1);
    }

    if (!ranges.empty())
      indices = parseRangeList<uint>(ranges);

  }

  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }

}



int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);
  
  parseOptions(argc, argv);

  AtomicGroup model = createSystem(model_name);

  AtomicGroup avg_subset = selectAtoms(model, avg_string);
  cerr << "Averaging over " << avg_subset.size() << " atoms.\n";

  pTraj traj = createTrajectory(traj_name, model);

  if (indices.empty()) {
    for (uint i = 0; i < traj->nframes(); ++i)
      indices.push_back(i);
  }

  cerr << "Using " << indices.size() << " frames from the trajectory...\n";

  // First, align...
  vector<XForm> xforms;
  if (align_string.empty()) {

    cerr << "Skipping alignment...\n";
    for (uint i=0; i<traj->nframes(); ++i)
      xforms.push_back(XForm());

  } else {

    AtomicGroup align_subset = selectAtoms(model, align_string);
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
  avgpdb.clearBonds();
  avgpdb.remarks().add(header);
  cout << avgpdb;
}
