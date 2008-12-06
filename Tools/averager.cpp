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


struct Globals {
  Globals() : trajmin(0), trajmax(0) { }

  string align_string;
  string avg_string;
  uint trajmin, trajmax;
  double alignment_tol;
  string model_name, traj_name;
};

Globals globals;



void parseOptions(int argc, char *argv[]) {

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("align,a", po::value<string>(&globals.align_string)->default_value("name == 'CA'"),"Align using this selection")
      ("average,A", po::value<string>(&globals.avg_string)->default_value("!(hydrogen || segid == 'SOLV' || segid == 'BULK')"), "Average over this selection")
      ("range,r", po::value<string>(), "Range of frames to average over (min:max)");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&globals.model_name), "Model filename")
      ("traj", po::value<string>(&globals.traj_name), "Trajectory filename");

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
      cerr << "Usage- rmsds [options] model-name trajectory-name\n";
      cerr << generic;
      exit(-1);
    }

    if (vm.count("range")) {
      string rangespec = vm["range"].as<string>();
      int i = sscanf(rangespec.c_str(), "%u:%u", &globals.trajmin, &globals.trajmax);
      if (i != 2) {
        cerr << "Error - invalid range specified for trajectory.\n";
        exit(-1);
      }

    }
  }

  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }

}


vector<XForm> doAlign(const AtomicGroup& subset, pTraj traj) {

  boost::tuple<vector<XForm>, greal, int> res = loos::iterativeAlignment(subset, traj, globals.alignment_tol, 100);
  vector<XForm> xforms = boost::get<0>(res);
  greal rmsd = boost::get<1>(res);
  int iters = boost::get<2>(res);

  cerr << "Subset alignment with " << subset.size()
       << " atoms converged to " << rmsd << " rmsd after "
       << iters << " iterations.\n";

  return(xforms);
}


int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);
  
  parseOptions(argc, argv);

  AtomicGroup model = loos::createSystem(globals.model_name);

  AtomicGroup align_subset = loos::selectAtoms(model, globals.align_string);
  cerr << "Aligning with " << align_subset.size() << " atoms.\n";

  AtomicGroup avg_subset = loos::selectAtoms(model, globals.avg_string);
  cerr << "Averaging over " << avg_subset.size() << " atoms.\n";

  pTraj traj = loos::createTrajectory(globals.traj_name, model);

  globals.trajmax = (globals.trajmax == 0) ? traj->nframes() : globals.trajmax+1;

  cerr << "Aligning...\n";
  vector<XForm> xforms = doAlign(align_subset, traj);
  cerr << "Averaging...\n";

  AtomicGroup avg = loos::averageStructure(avg_subset, xforms, traj);
  
  PDB avgpdb = PDB::fromAtomicGroup(avg);
  avgpdb.remarks().add(header);
  cout << avgpdb;
}
