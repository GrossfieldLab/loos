/*
  rmsds.cpp

  Pair-wise RMSD
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

typedef RealMatrix              Matrix;


const int matrix_precision = 2;    // Controls precision in output matrix

int verbosity;



// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() { }

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("noout,N", po::value<bool>(&noop)->default_value(false), "Do not output the matrix (i.e. only calc pair-wise RMSD stats)")
      ("sel1", po::value<string>(&sel1)->default_value("name == 'CA'"), "Atom selection for first system")
      ("skip1", po::value<uint>(&skip1)->default_value(0), "Skip n-frames of first trajectory")
      ("range1", po::value<string>(&range1), "Matlab-style range of frames to use from first trajectory")
      ("sel2", po::value<string>(&sel2)->default_value("name == 'CA'"), "Atom selection for second system")
      ("skip2", po::value<uint>(&skip2)->default_value(0), "Skip n-frames of second trajectory")
      ("range2", po::value<string>(&range2), "Matlab-style range of frames to use from second trajectory");

  }

  void addHidden(po::options_description& o) {
    o.add_options()
      ("model1", po::value<string>(&model1), "Model-1 Filename")
      ("traj1", po::value<string>(&traj1), "Traj-1 Filename")
      ("model2", po::value<string>(&model2), "Model-2 Filename")
      ("traj2", po::value<string>(&traj2), "Traj-2 Filename");
  }

  void addPositional(po::positional_options_description& pos) {
    pos.add("model1", 1);
    pos.add("traj1", 1);
    pos.add("model2", 1);
    pos.add("traj2", 1);
  }


  bool check(po::variables_map& m) {
    return( ! ( (m.count("model1") && m.count("traj1")) && !(m.count("model2") ^ m.count("traj2"))) );
  }


  string help() const {
    return("model-1 trajectory-1 [model-2 trajectory-2]");
  }


  string print() const {
    ostringstream oss;
    oss << boost::format("noout=%d,sel1='%s',skip1=%d,range1='%s',sel2='%s',skip2=%d,range2='%s',model1='%s',traj1='%s',model2='%s',traj2='%s'")
      % noop
      % sel1
      % skip1
      % range1
      % sel2
      % skip2
      % range2
      % model1
      % traj1
      % model2
      % traj2;

    return(oss.str());
  }


  bool noop;
  uint skip1, skip2;
  string range1, range2;
  string model1, traj1, model2, traj2;
  string sel1, sel2;
};


Matrix singleTrajectory(ToolOptions* topts) {
  AtomicGroup model = createSystem(topts->model1);
  pTraj traj = createTrajectory(topts->traj1, model);
  AtomicGroup subset = selectAtoms(model, topts->sel1);
  vector<uint> indices = assignTrajectoryFrames(traj, topts->range1, topts->skip1);

  PercentProgressWithTime watcher;
  PercentTrigger trigger(0.1);

  Matrix M(indices.size(), indices.size());
  double mean_rmsd = 0;
  double max_rmsd = 0;
  uint total = floor(indices.size()*indices.size()/2.0);

  ProgressCounter<PercentTrigger, EstimatingCounter> slayer(trigger, EstimatingCounter(total));
  if (verbosity > 0) {
    slayer.attach(&watcher);
    slayer.start();
  }

  AtomicGroup duplicate = subset.copy();
  for (uint j=1; j<indices.size(); ++j)
    for (uint i=0; i<j; ++i) {

      if (verbosity > 0)
        slayer.update();

      traj->readFrame(indices[j]);
      traj->updateGroupCoords(subset);

      traj->readFrame(indices[i]);
      traj->updateGroupCoords(duplicate);

      subset.alignOnto(duplicate);
      double r = subset.rmsd(duplicate);
      
      M(j, i) = r;
      M(i, j) = r;

      if (r > max_rmsd)
        max_rmsd = r;
      mean_rmsd += r;
    }

  if (verbosity > 0)
    slayer.finish();

  mean_rmsd /= total;
  cerr << boost::format("Max rmsd = %f, mean rmsd = %f\n") % max_rmsd % mean_rmsd;

  return(M);
}


Matrix twoTrajectories(ToolOptions* topts) {
  AtomicGroup model1 = createSystem(topts->model1);
  pTraj traj1 = createTrajectory(topts->traj1, model1);
  AtomicGroup subset1 = selectAtoms(model1, topts->sel1);
  vector<uint> indices1 = assignTrajectoryFrames(traj1, topts->range1, topts->skip1);

  AtomicGroup model2 = createSystem(topts->model2);
  pTraj traj2 = createTrajectory(topts->traj2, model2);
  AtomicGroup subset2 = selectAtoms(model2, topts->sel2);
  vector<uint> indices2 = assignTrajectoryFrames(traj2, topts->range2, topts->skip2);


  PercentProgressWithTime watcher;
  PercentTrigger trigger(0.1);

  Matrix M(indices1.size(), indices2.size());
  double mean_rmsd = 0;
  double max_rmsd = 0;
  uint total = floor(M.rows() * M.cols());

  ProgressCounter<PercentTrigger, EstimatingCounter> slayer(trigger, EstimatingCounter(total));
  if (verbosity > 0) {
    slayer.attach(&watcher);
    slayer.start();
  }

  for (uint j=0; j<indices1.size(); ++j)
    for (uint i=0; i<indices2.size(); ++i) {

      if (verbosity > 0)
        slayer.update();

      traj1->readFrame(indices1[j]);
      traj1->updateGroupCoords(model1);

      traj2->readFrame(indices2[i]);
      traj2->updateGroupCoords(model2);

      subset1.alignOnto(subset2);
      double r = subset1.rmsd(subset2);
      
      M(j, i) = r;
      if (r > max_rmsd)
        max_rmsd = r;
      mean_rmsd += r;
    }

  if (verbosity > 0)
    slayer.finish();

  mean_rmsd /= total;
  cerr << boost::format("Max rmsd = %f, mean rmsd = %f\n") % max_rmsd % mean_rmsd;

  return(M);
}


int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);
  
  opts::BasicOptions* bopts = new opts::BasicOptions;
  opts::BasicSelection* sopts = new opts::BasicSelection;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  verbosity = bopts->verbosity;
  Matrix M;

  if (topts->model2.empty() && topts->traj2.empty())
    M = singleTrajectory(topts);
  else
    M = twoTrajectories(topts);

  if (!topts->noop) {
    cout << "# " << header << endl;
    cout << setprecision(matrix_precision) << M;
  }
}


  

