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

using namespace std;
using namespace loos;


namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;




typedef Math::Matrix<double, Math::Triangular> Matrix;

// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() : tol(1e-6), iterate(false), noop(false) { }

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("iterative,I", po::value<bool>(&iterate)->default_value(iterate),"Use iterative alignment method")
      ("tolerance,T", po::value<double>(&tol)->default_value(tol), "Tolerance to use for iterative alignment")
      ("noout,N", po::value<bool>(&noop)->default_value(noop), "Do not output the matrix (i.e. only calc pair-wise RMSD stats)");
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("iterative=%d, tolerance=%f, noout=%d") % iterate % tol % noop;
    return(oss.str());
  }


  double tol;
  bool iterate;
  bool noop;
};


void doAlign(vector<AtomicGroup>& frames, const AtomicGroup& subset, pTraj traj, vector<uint>& idx, const double tol) {

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


Matrix interFrameRMSD(vector<AtomicGroup>& frames, bool iterate) {
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
  
  opts::BasicOptions* bopts = new opts::BasicOptions;
  opts::BasicSelection* sopts = new opts::BasicSelection;
  opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(tropts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  AtomicGroup molecule = tropts->model;
  pTraj ptraj = tropts->trajectory;
  AtomicGroup subset = selectAtoms(molecule, sopts->selection);
  cerr << "Selected " << subset.size() << " atoms in subset.\n";

  vector<uint> indices = tropts->frameList();

  vector<AtomicGroup> frames;
  if (topts->iterate) {
    cerr << "Aligning...\n";
    doAlign(frames, subset, ptraj, indices, topts->tol);
  } else
    readFrames(frames, subset, ptraj, indices);

  cerr << "Computing RMSD matrix...\n";
  Matrix M = interFrameRMSD(frames, topts->iterate);

  if (!topts->noop) {
    // Note:  using the operator<< on a matrix here will write it out as a full matrix
    //        i.e. not the special triangular format.
    cout << "# " << header << endl;
    cout << M;
  }
}


  

