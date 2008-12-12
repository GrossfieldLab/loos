/*
  rmsds.cpp

  Computes inter-frame rmsds.  Can run in two modes:
    o iterative
    o pair-wise

  Iterative computes an optimal global alignment through an interative
  scheme (see aligner.cpp for more details).  Pair-wise will compute
  the best pair-wise superposition prior to computing the RMSD.
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

namespace po = boost::program_options;
using namespace std;
using namespace loos;



typedef Math::Matrix<double, Math::Triangular> Matrix;


struct Globals {
  string model_name;
  string traj_name;
  string alignment;
  double tol;
  bool iterate;
};


Globals globals;

void parseOptions(int argc, char *argv[]) {

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("align,a", po::value<string>(&globals.alignment)->default_value("name == 'CA'"), "Align using this selection")
      ("iterative,i", po::value<bool>(&globals.iterate)->default_value(false),"Use iterative alignment method")
      ("tolerance,t", po::value<double>(&globals.tol)->default_value(1e-6), "Tolerance to use for iterative alignment");

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
  }

  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }

}


void doAlign(vector<AtomicGroup>& frames, const AtomicGroup& subset, pTraj traj) {
  uint n = traj->nframes();

  for (uint i = 0; i<n; i++) {
    AtomicGroup frame = subset.copy();
    traj->readFrame(i);
    traj->updateGroupCoords(frame);
    frames.push_back(frame);
  }

  boost::tuple<vector<XForm>, greal, int> res = iterativeAlignment(frames, globals.tol, 100);
  vector<XForm> xforms = boost::get<0>(res);
  greal rmsd = boost::get<1>(res);
  int iters = boost::get<2>(res);

  cerr << "Subset alignment with " << subset.size()
       << " atoms converged to " << rmsd << " rmsd after "
       << iters << " iterations.\n";
}


void readFrames(vector<AtomicGroup>& frames, const AtomicGroup& subset, pTraj traj) {
  uint n = traj->nframes();

  for (uint i=0; i<n; i++) {
    AtomicGroup frame = subset.copy();
    traj->readFrame(i);
    traj->updateGroupCoords(frame);
    frames.push_back(frame);
  }
}


Matrix interFrameRMSD(vector<AtomicGroup>& frames) {
  uint n = frames.size();
  Matrix M(n, n);
  uint i;

  uint total = n*(n+1)/2;
  uint delta = total / 4;
  uint k = 0;
  uint j;

  for (j=0; j<n; j++) {
    AtomicGroup jframe = frames[j].copy();
    for (i=0; i<=j; i++, k++) {
      double rmsd;

      if (globals.iterate)
        rmsd = jframe.rmsd(frames[i]);
      else {
        (void)jframe.alignOnto(frames[i]);
        rmsd = jframe.rmsd(frames[i]);
      }

      M(j,i) = rmsd;
      
      if (k % delta == 0) {
        float percent = k * 100.0 / total;
        cerr << setprecision(3) << percent << "% complete\n";
      }
    }
  }

  return(M);
}


int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  AtomicGroup molecule = createSystem(globals.model_name);
  pTraj ptraj = createTrajectory(globals.traj_name, molecule);
  AtomicGroup subset = selectAtoms(molecule, globals.alignment);
  cerr << "Selected " << subset.size() << " atoms.\n";

  vector<AtomicGroup> frames;
  if (globals.iterate) {
    cerr << "Aligning...\n";
    doAlign(frames, subset, ptraj);
  } else
    readFrames(frames, subset, ptraj);

  Matrix M = interFrameRMSD(frames);

  // Note:  using the operator<< on a matrix here will write it out as a full matrix
  //        i.e. not the special triangular format.
  cout << "# " << header << endl;
  cout << M;
}


  

