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
#include <boost/format.hpp>

namespace po = boost::program_options;
using namespace std;
using namespace loos;



typedef Math::Matrix<double, Math::Triangular> Matrix;

const int matrix_precision = 2;    // Controls precision in output matrix



string model_name;
string traj_name;
string alignment;
double tol;
bool iterate;
bool no_output;
vector<uint> indices;


void parseOptions(int argc, char *argv[]) {

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("align,a", po::value<string>(&alignment)->default_value("name == 'CA'"), "Align using this selection")
      ("iterative,i", po::value<bool>(&iterate)->default_value(false),"Use iterative alignment method")
      ("tolerance,t", po::value<double>(&tol)->default_value(1e-6), "Tolerance to use for iterative alignment")
      ("range,r", po::value< vector<string> >(), "Frames of the trajectory to use (list of Octave-style ranges)")
      ("noout,n", po::value<bool>(&no_output)->default_value(false), "Do not output the matrix (i.e. only calc pair-wise RMSD stats)");

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
      cerr << "Usage- rmsds [options] model-name trajectory-name\n";
      cerr << generic;
      exit(-1);
    }
    
    if (vm.count("range")) {
      vector<string> ranges = vm["range"].as< vector<string> >();
      indices = parseRangeList<uint>(ranges);
    }


  }

  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }

}


void doAlign(vector<AtomicGroup>& frames, const AtomicGroup& subset, pTraj traj, vector<uint>& idx) {

  for (vector<uint>::iterator i = idx.begin(); i != idx.end(); ++i) {
    AtomicGroup frame = subset.copy();
    traj->readFrame(*i);
    traj->updateGroupCoords(frame);
    frames.push_back(frame);
  }

  boost::tuple<vector<XForm>, greal, int> res = iterativeAlignment(frames, tol, 100);
  vector<XForm> xforms = boost::get<0>(res);
  greal rmsd = boost::get<1>(res);
  int iters = boost::get<2>(res);

  cerr << "Subset alignment with " << subset.size()
       << " atoms converged to " << rmsd << " rmsd after "
       << iters << " iterations.\n";
}


void readFrames(vector<AtomicGroup>& frames, const AtomicGroup& subset, pTraj traj, vector<uint>& idx) {

  for (vector<uint>::iterator i = idx.begin(); i != idx.end(); ++i) {
    AtomicGroup frame = subset.copy();
    traj->readFrame(*i);
    traj->updateGroupCoords(frame);
    frames.push_back(frame);
  }
}


Matrix interFrameRMSD(vector<AtomicGroup>& frames) {
  uint n = frames.size();
  Matrix M(n, n);
  uint i;

  uint k = 0;
  uint j;

  double max = 0.0;
  double mean = 0.0;
  uint total = n*(n+1)/2;

  PercentProgressWithTime watcher;
  PercentTrigger trigger(0.25);

  ProgressCounter<PercentTrigger, EstimatingCounter> slayer(trigger, EstimatingCounter(total));
  slayer.attach(&watcher);
  slayer.start();

  for (j=0; j<n; j++) {
    AtomicGroup jframe = frames[j].copy();
    for (i=0; i<=j; i++, k++) {
      double rmsd;

      slayer.update();

      if (iterate)
        rmsd = jframe.rmsd(frames[i]);
      else {
        (void)jframe.alignOnto(frames[i]);
        rmsd = jframe.rmsd(frames[i]);
      }

      M(j,i) = rmsd;

      if (j != i) {
        mean += rmsd;
        if (rmsd > max)
          max = rmsd;
      }
      
    }
  }

  slayer.finish();
  mean /= total;
  cerr << boost::format("Max rmsd = %f, mean rmsd = %f\n") % max % mean;

  return(M);
}


int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  AtomicGroup molecule = createSystem(model_name);
  pTraj ptraj = createTrajectory(traj_name, molecule);
  AtomicGroup subset = selectAtoms(molecule, alignment);
  cerr << "Selected " << subset.size() << " atoms in subset.\n";

  if (indices.empty())
    for (uint i=0; i<ptraj->nframes(); ++i)
      indices.push_back(i);

  vector<AtomicGroup> frames;
  if (iterate) {
    cerr << "Aligning...\n";
    doAlign(frames, subset, ptraj, indices);
  } else
    readFrames(frames, subset, ptraj, indices);

  cerr << "Computing RMSD matrix...\n";
  Matrix M = interFrameRMSD(frames);

  if (!no_output) {
    // Note:  using the operator<< on a matrix here will write it out as a full matrix
    //        i.e. not the special triangular format.
    cout << "# " << header << endl;
    cout << setprecision(matrix_precision) << M;
  }
}


  

