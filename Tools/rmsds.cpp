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
#include <unistd.h>

using namespace std;
using namespace loos;


namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

typedef RealMatrix              Matrix;

#if defined(__linux__) || defined(__APPLE__)
const int matrix_precision = 2;    // Controls precision in output matrix
#endif

int verbosity;


// If the estimated cache memory is more than this fraction of physical memory,
// issue a warning to the user to consider turning off the cache
const double cache_memory_fraction_warning = 0.66;



// @cond TOOLS_INTERNAL

string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\n"
    "\tCalculate a pair-wise RMSD for a trajectory (or two trajectories)\n"
    "DESCRIPTION\n"
    "\n"
    "\tThis tool calculates the pair-wise RMSD between each structure in a trajectory\n"
    "or, alternatively, between each structure in two different trajectories.  In the single\n"
    "trajectory case, the ith structure is aligned with the jth structure and the RMSD calculated.\n"
    "This is stored in a matrix, i.e. R(j, i) = d(S_i, S_j).  The block-structure is indicative\n"
    "of sets of similar conformations.  The presence (or lack thereof) of multiple cross-peaks\n"
    "is diagnostic of the sampling quality of a simulation.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\trmsds model.pdb simulation.dcd >rmsd.asc\n"
    "This example uses all alpha-carbons and every frame in the trajectory.\n"
    "\n"
    "\trmsds inactive.pdb inactive.dcd active.pdb active.dcd >rmsd.asc\n"
    "This example uses all alpha-carbons and compares the \"inactive\" simulation\n"
    "with the \"active\" one.\n"
    "\n"
    "\trmsds --sel1 'resid <= 100 && name == \"CA\"' model.pdb simulation.dcd >rmsds.asc\n"
    "This example calculates the pair-wise RMSD using only the first 100 alpha-carbons\n"
    "\n"
    "\trmsds --sel1 'resid <= 50 && name == \"CA\"' \\\n"
    "\t  --sel2 'resid >=20 && resid <= 69 && name == \"CA\"' \\\n"
    "\t  inactive.pdb inactive.dcd active.pdb active.dcd >rmsd.asc\n"
    "This example compares two trajectories, active and inactive, and uses different selections\n"
    "for both: the first 50 residues from the inactive and residues 20-69 from the active.\n"
    "\n"
    "NOTES\n"
    "\tWhen using two trajectories, the selections must match both in number of atoms selected\n"
    "and in the sequence of atoms (i.e. the first atom in the --sel2 selection is\n" 
    "matched with the first atom in the --sel2 selection.)\n"
    "\n"
    "SEE ALSO\n"
    "\trmsd2ref\n"
    "\n";

  return(msg);
}



class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() { }

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("cache,C", po::value<bool>(&cache)->default_value(true), "Cache frames from the trajectory for speed")
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
    oss << boost::format("cached=%d,noout=%d,sel1='%s',skip1=%d,range1='%s',sel2='%s',skip2=%d,range2='%s',model1='%s',traj1='%s',model2='%s',traj2='%s'")
      % cache
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


  bool cache;
  bool noop;
  uint skip1, skip2;
  string range1, range2;
  string model1, traj1, model2, traj2;
  string sel1, sel2;
};



// This is the top-level interface for getting frames from the trajectory.
// It allows the user to select whether the frames will be cached in memory
// or always read from disk

class TrajInterface 
{
public:
  TrajInterface(const AtomicGroup& model, pTraj& traj) : _model(model.copy()), _traj(traj) 
  {
  }
  

  virtual ~TrajInterface() 
  {
  }
  

  // Return the ith structure
  virtual AtomicGroup getStructure(const uint i) =0;

  // Return the total number of structures
  virtual uint size() const =0;

protected:
  AtomicGroup _model;
  pTraj _traj;
};


class CachedTraj : public TrajInterface 
{
public:
  CachedTraj(const AtomicGroup& model, pTraj& traj, const vector<uint>& frames) :
    TrajInterface(model, traj)
  {
    _model.clearBonds();    // Save a little space

#if defined(__linux__) || defined(__APPLE__)
    // Estimate spaced consumed by cache and warn if it's large (and potentially swapping)
    long ms = estimateModelSize(_model);
    long ts = estimateEnsembleSize(frames.size());
    long mem = availMemory();
    double used = static_cast<double>(ms * ts) / mem;
    
    if (used >= cache_memory_fraction_warning) {
      long mb = (ms * ts) >> 20;

      cerr << boost::format("***WARNING***\nThe estimated trajectory cache size is %.1f%% of your memory (%d MB)\n")
	% (used * 100.0)
	% mb;

      cerr << "If your machine starts swapping, try using the --cache=0 option on the command line\n";
    }
#endif
    
    readTrajectory(_ensemble, model, traj);
    
  }


#if defined(__linux__) || defined(__APPLE__)

  // Should consider using _SC_AVPHYS_PAGES instead?
  long availMemory() const
  {
    long pagesize = sysconf(_SC_PAGESIZE);
    long pages = sysconf(_SC_PHYS_PAGES);

    return(pagesize * pages);
  }
  

  // Manually count size of a model (including contained strings)
  long estimateModelSize(const AtomicGroup& model) const
  {
    long s = sizeof(model);
    for (uint i=0; i<model.size(); ++i) {
      s += sizeof(Atom);
      s += model[i]->name().size();
      s += model[i]->altLoc().size();
      s += model[i]->chainId().size();
      s += model[i]->resname().size();
      s += model[i]->segid().size();
      s += model[i]->iCode().size();
      s += model[i]->PDBelement().size();
      s += model[i]->recordName().size();
      s += 8;
    }
    
    return(s);
    }
    

  // Assume vector capacity will always be a power of 2
  long estimateEnsembleSize(const uint n) const
  {
    long actual = 2;

    while (actual < n)
      actual <<= 1;
    
    return(actual);
  }
  
#endif
  
  // Each AtomicGroup in the vector is a copy, so we can get away with just
  // returning the element itself
  AtomicGroup getStructure(const uint i) 
  {
    return(_ensemble[i]);
  }
  
  uint size() const 
  {
    return(_ensemble.size());
  }
  
private:
  vector<AtomicGroup> _ensemble;
};



class NonCachedTraj : public TrajInterface 
{
public:
  NonCachedTraj(const AtomicGroup& model, pTraj& traj, const vector<uint>& frames) :
    TrajInterface(model, traj),
    _frames(frames)
  {
  }
  

  // Must return a copy here so that the atoms are not shared...
  AtomicGroup getStructure(const uint i) 
  {
    AtomicGroup frame = _model.copy();
    
    _traj->readFrame(_frames[i]);
    _traj->updateGroupCoords(frame);
    return(frame);
  }
  
  uint size() const 
  {
    return(_frames.size());
  }
  
private:
  vector<uint> _frames;
};




// @endcond TOOLS_INTERNAL


Matrix singleTrajectory(ToolOptions* topts) {
  AtomicGroup model = createSystem(topts->model1);
  pTraj traj = createTrajectory(topts->traj1, model);
  AtomicGroup subset = selectAtoms(model, topts->sel1);
  vector<uint> indices = assignTrajectoryFrames(traj, topts->range1, topts->skip1);

  TrajInterface *traji;
  if (topts->cache)
    traji = new CachedTraj(subset, traj, indices);
  else
    traji = new NonCachedTraj(subset, traj, indices);
  

  PercentProgressWithTime watcher;
  PercentTrigger trigger(0.1);

  Matrix M(traji->size(), traji->size());
  double mean_rmsd = 0;
  double max_rmsd = 0;
  uint total = floor(traji->size()*traji->size()/2.0);

  ProgressCounter<PercentTrigger, EstimatingCounter> slayer(trigger, EstimatingCounter(total));
  if (verbosity > 0) {
    slayer.attach(&watcher);
    slayer.start();
  }

  AtomicGroup duplicate = subset.copy();
  for (uint j=1; j<traji->size(); ++j) {

    AtomicGroup model_j = traji->getStructure(j);

    for (uint i=0; i<j; ++i) {

      if (verbosity > 0)
        slayer.update();


      AtomicGroup model_i = traji->getStructure(i);

      model_i.alignOnto(model_j);
      double r = model_j.rmsd(model_i);
      
      M(j, i) = r;
      M(i, j) = r;

      if (r > max_rmsd)
        max_rmsd = r;
      mean_rmsd += r;
    }
  }
  

  if (verbosity > 0)
    slayer.finish();

  mean_rmsd /= total;
  cerr << boost::format("Max rmsd = %f, mean rmsd = %f\n") % max_rmsd % mean_rmsd;

  delete traji;

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
  
  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(topts);
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


  

