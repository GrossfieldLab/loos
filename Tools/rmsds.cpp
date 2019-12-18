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
#include <boost/thread/thread.hpp>


using namespace std;
using namespace loos;


typedef boost::thread*   Thread;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;



// const int matrix_precision = 2;    // Controls precision in output matrix

int verbosity;


// If the estimated cache memory is more than this fraction of physical memory,
// issue a warning to the user to consider turning off the cache
// Note: the total app size may be 20-30% larger than the cache estimate, so
//       take that into consideration when setting the warning threshold

const double cache_memory_fraction_warning = 0.66;

long used_memory = 0;
bool report_stats = false;



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
    "\tThe requested subset for each frame is cached in memory for better performance.\n"
    "If the memory used by the cache gets too large, your machine may swap and dramatically slow\n"
    "down.  The tool will try to warn you if this is a possibility.  To use less memory, subsample\n"
    "the trajectory either by using the --range1 and --range2 options, or use subsetter to pre-process\n"
    "the trajectory.\n"
    "\n"
    "\tThis tool can be run in parallel with multiple threads for performance.  The --threads option\n"
    "controls how many threads are used.  The default is 1 (non-parallel).  Setting it to 0 will use\n"
    "as many threads as possible.  Note that if LOOS was built using a multi-threaded math library,\n"
    "then some care should be taken in how many threads are used for this tool, though it is unlikely\n"
    "that there will be a conflict.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\trmsds model.pdb simulation.dcd >rmsd.asc\n"
    "This example uses all alpha-carbons and every frame in the trajectory.\n"
    "\n"
    "\trmsds --threads=8 model.pdb simulation.dcd >rmsd.asc\n"
    "This example uses all alpha-carbons and every frame in the trajectory, run\n"
    "in parallel with 8 threads of execution.\n"
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
      ("noout,N", po::value<bool>(&noop)->default_value(false), "Do not output the matrix (i.e. only calc pair-wise RMSD stats)")
      ("threads", po::value<uint>(&nthreads)->default_value(1), "Number of threads to use (0=all available)")
      ("sel1", po::value<string>(&sel1)->default_value("name == 'CA'"), "Atom selection for first system")
      ("skip1", po::value<uint>(&skip1)->default_value(0), "Skip n-frames of first trajectory")
      ("range1", po::value<string>(&range1), "Matlab-style range of frames to use from first trajectory")
      ("sel2", po::value<string>(&sel2)->default_value("name == 'CA'"), "Atom selection for second system")
      ("skip2", po::value<uint>(&skip2)->default_value(0), "Skip n-frames of second trajectory")
      ("range2", po::value<string>(&range2), "Matlab-style range of frames to use from second trajectory")
      ("stats", po::value<bool>(&stats)->default_value(false), "Show some statistics for matrix")
      ("precision,p", po::value<uint>(&matrix_precision)->default_value(2), "Write out matrix coefficients with this many digits.");
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
    oss << boost::format("stats=%d,matrix_precision=%d,noout=%d,nthreads=%d,sel1='%s',skip1=%d,range1='%s',sel2='%s',skip2=%d,range2='%s',model1='%s',traj1='%s',model2='%s',traj2='%s'")
      % stats
      % matrix_precision
      % noop
      % nthreads
      % matrix_precision
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


  bool stats;
  bool noop;
  uint skip1, skip2;
  uint nthreads;
  uint matrix_precision;
  string range1, range2;
  string model1, traj1, model2, traj2;
  string sel1, sel2;
};

typedef vector<double>    vecDouble;
typedef vector<vecDouble>   vMatrix;


// @endcond TOOLS_INTERNAL


// --------------------------------------------------------------------------------------

// Parcels out work to the compute threads...  Work is given to the threads
// one row at a time.

class Master {
public:

  Master(const uint nr, const bool tr, const bool b) : _toprow(0), _maxrow(nr), _updatefreq(500), _triangle(tr),
                                                       _verbose(b), _start_time(time(0))
  {
    if (_triangle)
      _total = _maxrow*(_maxrow-1) / 2;
    else
      _total = _maxrow;
  }

  // Checks whether there are any columns left to work on
  // and places the column index into the passed pointer.

  bool workAvailable(uint* ip) 
  {

    _mtx.lock();
    if (_toprow >= _maxrow) {
      _mtx.unlock();
      return(false);
    }
    *ip = _toprow++;

    if (_verbose)
      if (_toprow % _updatefreq == 0)
        updateStatus();
    
    _mtx.unlock();
    return(true);
  }


  void updateStatus() {
    time_t dt = elapsedTime();
    uint work_done = _triangle ? (_toprow * (_toprow-1) / 2) : (_toprow);
    uint work_left = _total - work_done;
    uint d = work_left * dt / work_done;    // rate = work_done / dt;  d = work_left / rate;
    
    uint hrs = d / 3600;
    uint remain = d % 3600;
    uint mins = remain / 60;
    uint secs = remain % 60;
    
    cerr << boost::format("Row %5d /%5d, Elapsed = %5d s, Remaining = %02d:%02d:%02d\n")
      % _toprow % _maxrow % dt % hrs % mins % secs;
  }


  
  time_t elapsedTime() const 
  {
    return (time(0) - _start_time);
  }


private:
  uint _toprow, _maxrow;
  uint _updatefreq;
  bool _triangle;
  bool _verbose;
  time_t _start_time;
  uint _total;
  boost::mutex _mtx;

};






/*
  Worker thread processes a column of the all-to-all matrix.  Gets which
  column to work on from the associated Master object.
*/


// Worker for two different trajectories

class DualWorker 
{
public:
  DualWorker(RealMatrix* R, vMatrix* T1, vMatrix* T2, Master* M) : _R(R), _T1(T1), _T2(T2), _M(M), _maxcol(T2->size()) { }


  DualWorker(const DualWorker& w) 
  {
    _R = w._R;
    _T1 = w._T1;
    _T2 = w._T2;
    _M = w._M;
    _maxcol = w._maxcol;
  }
  

  void calc(const uint i) 
  {
    for (uint j=0; j<_maxcol; ++j) {
      double d = loos::alignment::centeredRMSD((*_T1)[i], (*_T2)[j]);
      (*_R)(i, j) = d;
    }
  }

  void operator()() 
  {
    uint i;
    
    while (_M->workAvailable(&i))
      calc(i);
  }
  

private:
  RealMatrix* _R;
  vMatrix* _T1;
  vMatrix* _T2;
  Master* _M;
  uint _maxcol;
};



// Worker for self all-to-all

class SingleWorker 
{
public:
  SingleWorker(RealMatrix* R, vMatrix* T, Master* M) : _R(R), _T(T), _M(M) { }


  SingleWorker(const SingleWorker& w) 
  {
    _R = w._R;
    _T = w._T;
    _M = w._M;
  }
  

  void calc(const uint i) 
  {
    for (uint j=0; j<i; ++j) {
      double d = loos::alignment::centeredRMSD((*_T)[i], (*_T)[j]);
      (*_R)(j, i) = (*_R)(i, j) = d;
    }
  }

  void operator()() 
  {
    uint i;
    
    while (_M->workAvailable(&i))
      calc(i);
  }
  

private:
  RealMatrix* _R;
  vMatrix* _T;
  Master* _M;
};



// Top level object/interface.  Will create np Worker threads, cloned from the
// one passed into the constructor.

template<class W>
class Threader 
{
public:
  Threader(W* wrkr, const uint np) :
    worker(wrkr),
    threads(vector<Thread>(np))
  {

    for (uint i=0; i<np; ++i) {
      W w(*worker);
      threads[i] = new boost::thread(w);
    }
    
  }


  void join() 
  {
    for (uint i=0; i<threads.size(); ++i)
      threads[i]->join();
  }
  

  ~Threader() 
  {
    for (uint i=0; i<threads.size(); ++i)
      delete threads[i];
  }

  W* worker;
  vector<Thread> threads;
};


// --------------------------------------------------------------------------------------


void showStatsHalf(const RealMatrix& R) {
  uint total = (R.rows() * (R.rows()-1)) / 2; 

  double avg = 0.0;
  double max = 0.0;
  for (uint j=1; j<R.rows(); ++j)
    for (uint i=0; i<j; ++i) {
      avg += R(j, i);
      if (R(j, i) > max)
        max = R(j, i);
    }
    
  avg /= total;
  cerr << boost::format("Max rmsd = %.4f, avg rmsd = %.4f\n") % max % avg;
}


void showStatsWhole(const RealMatrix& R) {
  uint total = R.rows() * R.cols();

  double avg = 0.0;
  double max = 0.0;
  for (ulong i=0; i<total; ++i) {
    avg += R[i];
    if (R[i] > max)
      max = R[i];
  }
    
  avg /= total;

  cerr << boost::format("Max rmsd = %.4f, avg rmsd = %.4f\n") % max % avg;
}



void centerTrajectory(alignment::vecMatrix& U) {
  for (uint i=0; i<U.size(); ++i)
    alignment::centerAtOrigin(U[i]);
}



void checkMemoryUsage(long mem) {
  if (!mem)
    return;

  double used = static_cast<double>(used_memory) / mem;

  if (verbosity > 2)
    cerr << boost::format("Memory: available=%d GB, estimated used=%.2f MB\n")
      % (mem >> 30) % (static_cast<double>(used_memory) / (1lu<<20) );
    
  if (used >= cache_memory_fraction_warning) {
    cerr << boost::format("***WARNING***\nThe estimated memory used is %.1f%% (%d MB) of your total memory (%d GB).\n")
      % (used * 100.0)
      % (used_memory >> 20)
      % (mem >> 30);
    
    cerr << "If your machine starts swapping, try subsampling the trajectories\n";
  }
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
  report_stats = (verbosity || topts->noop);
  AtomicGroup model = createSystem(topts->model1);
  pTraj traj = createTrajectory(topts->traj1, model);
  AtomicGroup subset = selectAtoms(model, topts->sel1);
  vector<uint> indices = assignTrajectoryFrames(traj, topts->range1, topts->skip1);

  long mem = availableMemory();
  uint nthreads = topts->nthreads ? topts->nthreads : boost::thread::hardware_concurrency();
  
  if (verbosity > 1) {
    cerr << "Using " << nthreads << " threads\n";
    cerr << "Reading trajectory - " << topts->traj1 << endl;
  }
  vMatrix T = readCoords(subset, traj, indices, verbosity > 1);
  used_memory += T.size() * T[0].size() * sizeof(vMatrix::value_type::value_type);   // Coords matrix
  used_memory += T.size() * T.size() * sizeof(RealMatrix::element_type);             // RMSDS matrix
  checkMemoryUsage(mem);
  centerTrajectory(T);

  RealMatrix M;
  if (topts->model2.empty()) {

    if (verbosity > 1)
      cerr << "Calculating RMSD...\n";
    M = RealMatrix(T.size(), T.size());
    Master master(T.size(), true, verbosity);
    SingleWorker worker(&M, &T, &master);
    Threader<SingleWorker> threads(&worker, nthreads);
    threads.join();
    if (verbosity) 
      master.updateStatus();
    
    if (verbosity || topts->noop || topts->stats)
      showStatsHalf(M);
    
  } else {
    AtomicGroup model2 = createSystem(topts->model2);
    pTraj traj2 = createTrajectory(topts->traj2, model2);
    AtomicGroup subset2 = selectAtoms(model2, topts->sel2);
    vector<uint> indices2 = assignTrajectoryFrames(traj2, topts->range2, topts->skip2);

    if (verbosity > 1)
      cerr << "Reading trajectory - " << topts->traj2 << endl;
    vMatrix T2 = readCoords(subset2, traj2, indices2, verbosity > 1);
    used_memory += T2.size() * T2[0].size() * sizeof(double);
    checkMemoryUsage(mem);
    centerTrajectory(T2);

    if (verbosity > 1)
      cerr << "Calculating RMSD...\n";
    M = RealMatrix(T.size(), T2.size());
    Master master(T.size(), false, verbosity);
    DualWorker worker(&M, &T, &T2, &master);
    Threader<DualWorker> threads(&worker, nthreads);
    threads.join();

    if (verbosity)
      master.updateStatus();

    if (verbosity || topts->noop || topts->stats)
      showStatsWhole(M);
  }

  if (!topts->noop) {
    cout << "# " << header << endl;
    cout << setprecision(topts->matrix_precision) << M;
  }

}


  

