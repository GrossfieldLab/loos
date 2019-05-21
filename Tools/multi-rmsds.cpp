/*
  multi-rmsds.cpp

  RMSDS between multiple trajectories
*/



/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008,2016 Tod D. Romo
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
    "\tCalculate a pair-wise RMSD for multiple trajectories\n"
    "DESCRIPTION\n"
    "\n"
    "\tThis tool calculates the pair-wise RMSD between each structure in a multi-trajectory\n"
    "The ith structure is aligned with the jth structure and the RMSD calculated.\n"
    "This is stored in a matrix, i.e. R(j, i) = d(S_i, S_j).  The block-structure is indicative\n"
    "of sets of similar conformations.  The presence (or lack thereof) of multiple cross-peaks\n"
    "is diagnostic of the sampling quality of a simulation.  Cross-peaks between sub-blocks indicates\n"
    "similar conformations in multiple trajectories.\n"
    "\n"
    "\tThe requested subset for each frame is cached in memory for better performance.\n"
    "If the memory used by the cache gets too large, your machine may swap and dramatically slow\n"
    "down.  The tool will try to warn you if this is a possibility.  To use less memory, subsample\n"
    "the trajectory by using --skip or --stride, or use subsetter to pre-process the trajectory.\n"
    "\n"
    "\tThis tool can be run in parallel with multiple threads for performance.  The --threads option\n"
    "controls how many threads are used.  The default is 1 (non-parallel).  Setting it to 0 will use\n"
    "as many threads as possible.  Note that if LOOS was built using a multi-threaded math library,\n"
    "then some care should be taken in how many threads are used for this tool, though it is unlikely\n"
    "that there will be a conflict.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\tmulti-rmsds model.pdb sim1.dcd sim2.dcd sim3.dcd >rmsd.asc\n"
    "This example uses all alpha-carbons and every frame from each trajectory.\n"
    "\n"
    "\tmulti-rmsds --threads=8 model.pdb sim1.dcd sim2.dcd sim3.dcd >rmsd.asc\n"
    "This example uses all alpha-carbons and every frame in the trajectories, run\n"
    "in parallel with 8 threads of execution.\n"
    "\n"
    "\tmulti-rmsds --selection backbone --skip=50 --stride=10 model.pdb sim1.dcd sim2.dcd sim3.dcd >rmsds.asc\n"
    "This example uses the backbone atoms, and skips the first 50 frames from each trajectory,\n"
    "and only takes every 10th subsequent frame from each trajectory.\n"
    "\n"
    "SEE ALSO\n"
    "\trmsds, rmsd2ref, rms-overlap\n"
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
      ("stats", po::value<bool>(&stats)->default_value(false), "Show some statistics for matrix")
      ("precision,p", po::value<uint>(&matrix_precision)->default_value(2), "Write out matrix coefficients with this many digits.");
  }



  string print() const {
    ostringstream oss;
    oss << boost::format("stats=%d,noout=%d,nthreads=%d,matrix_precision=%d")
      % stats
      % noop
      % nthreads
      % matrix_precision;

    return(oss.str());
  }


  bool stats;
  bool noop;
  uint nthreads;
  uint matrix_precision;
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
  opts::BasicSelection* sopts = new opts::BasicSelection("name == 'CA'");
  opts::MultiTrajOptions* mtopts = new opts::MultiTrajOptions;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(mtopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  verbosity = bopts->verbosity;
  report_stats = (verbosity || topts->noop);
  AtomicGroup model = mtopts->model;
  pTraj traj = mtopts->trajectory;
  AtomicGroup subset = selectAtoms(model, sopts->selection);
  vector<uint> indices = mtopts->frameList();

  long mem = availableMemory();
  uint nthreads = topts->nthreads ? topts->nthreads : boost::thread::hardware_concurrency();

  if (verbosity > 1)
    cerr << "Using " << nthreads << " threads\n";

  vMatrix T = readCoords(subset, traj, indices, verbosity > 1);
  used_memory += T.size() * T[0].size() * sizeof(vMatrix::value_type::value_type);   // Coords matrix
  used_memory += T.size() * T.size() * sizeof(RealMatrix::element_type);             // RMSDS matrix
  checkMemoryUsage(mem);
  centerTrajectory(T);

  RealMatrix M;
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


  if (!topts->noop) {
    cout << "# " << header << endl;
    cout << mtopts->trajectoryTable();
    cout << setprecision(topts->matrix_precision) << M;
  }

}
