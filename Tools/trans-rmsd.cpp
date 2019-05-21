/*
  trans-rmsd.cpp

  RMSDS between two sets of multiple trajectories
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
#include <boost/tuple/tuple.hpp>
#include <boost/thread/thread.hpp>
#include <boost/algorithm/string.hpp>


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
    "XXX";

  return(msg);
}

class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() { }

  void addGeneric(po::options_description& o) {
    std::string modeltypes = "Model types:\n" + availableSystemFileTypes();
    std::string trajtypes = "Trajectory types:\n" + availableTrajectoryFileTypes();
    o.add_options()
      ("modeltype", po::value<std::string>(), modeltypes.c_str())
      ("set-A,A", po::value<string>(&set_A), "Space separated set of trajectories to compare pair-wise to B")
      ("set-B,B", po::value<string>(&set_B), "Space separated set of trajectories to compare pair-wise to A")
      ("skip,k", po::value<uint>(&skip)->default_value(0), "Number of frames to skip in sub-trajectories")
      ("stride,i", po::value<uint>(&stride)->default_value(1), "Step through sub-trajectories by this amount")
      ("range,r", po::value<std::string>(&frame_index_spec), "Which frames to use in composite trajectory")
      ("noout,N", po::value<bool>(&noop)->default_value(false), "Do not output the matrix (i.e. only calc pair-wise RMSD stats)")
      ("threads", po::value<uint>(&nthreads)->default_value(1), "Number of threads to use (0=all available)")
      ("cutoff,c", po::value<float>(&cutoff)->default_value(-1.0), "Outputs fraction of frame-pairs below cutoff.")
      ("stats", po::value<bool>(&stats)->default_value(false), "Show some statistics for matrix")
      ("precision,p", po::value<uint>(&matrix_precision)->default_value(2), "Write out matrix coefficients with this many digits.");
  }
  void addHidden(po::options_description& opts) {
    opts.add_options()
      ("model", po::value<std::string>(&model_name), "Model filename");
      //("traj", po::value< std::vector< std::string > >(&traj_names), "Trajectory filenames");
  }

  void addPositional(po::positional_options_description& pos) {
    pos.add("model", 1);
    //pos.add("traj_names", -1)
  }
  bool check(po::variables_map& map) {
    return( model_name.empty() || set_A.empty() || set_B.empty());
  }


  bool postConditions(po::variables_map& map) {
    if (map.count("modeltype")) {
      model_type = map["modeltype"].as<std::string>();
      model = createSystem(model_name, model_type);
    } else
      model = createSystem(model_name);
    // boost::split(results_vector, input_string, function_to_split_at(char c))
    boost::split(trajlist_A, set_A, boost::is_any_of(" "));
    boost::split(trajlist_B, set_B, boost::is_any_of(" "));
    mtrajA = MultiTrajectory(trajlist_A, model, skip, stride);
    mtrajB = MultiTrajectory(trajlist_B, model, skip, stride);
    trajectory_A = pTraj(&mtrajA, boost::lambda::_1);
    trajectory_B = pTraj(&mtrajB, boost::lambda::_1);

    return true;
  }


  std::vector<uint> frameList(pTraj trajectory) const {
    std::vector<uint> indices = assignTrajectoryFrames(trajectory, frame_index_spec, 0, 1);
    return uniquifyVector(indices);
  }

  // There was an analogous call in multitrajectoryoptions, leaving it out for now to get other things working.
  //std::string help() const { return("model trajectory [trajectory ...]"); }
  
  std::string trajectoryTable(MultiTrajectory mtraj) const {
    std::ostringstream oss;
    if (!frame_index_spec.empty())
      oss << "# Note- composite frame range used was '" << frame_index_spec << "'\n";
    oss << "# traj\tstart\tend\tfilename\n";
    uint start_cnt = 0;
    uint j = 0;
    for (uint i=0; i<mtraj.size(); ++i) {
      uint n = mtraj.nframes(i);
      if (n == 0)
        oss << boost::format("# Warning- '%s' was skipped due to insufficient frames\n") % mtraj[i]->filename();
      else {
        oss << boost::format("# %d\t%d\t%d\t%s\n")
          % j
          % start_cnt
          % (start_cnt + n - 1)
          % mtraj[i]->filename();
        ++j;
      }
      start_cnt += n;
    }

    return oss.str();
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("model='%s', modeltype='%s', skip=%d, stride=%d, trajlist_A=(")
      % model_name % model_type % skip % stride;
    for (uint i=0; i<trajlist_A.size(); i++)
      oss << "'" << trajlist_A[i] << "'" << (i < trajlist_A.size() -1 ? "," : "");
    oss << "), trajlist_B=(";
    for (uint i=0; i<trajlist_B.size(); ++i)
      oss << "'" << trajlist_B[i] << "'" << (i < trajlist_B.size() -1 ? "," : "");
    oss << ")";
    oss << boost::format("stats=%d,noout=%d,nthreads=%d,matrix_precision=%d")
      % stats
      % noop
      % nthreads
      % matrix_precision; 
    return(oss.str());
  }

  bool stats;
  bool noop;
  float cutoff;
  uint nthreads;
  std::string set_A;
  std::string set_B;
  std::vector<string> trajlist_A, trajlist_B;
  uint skip;
  uint stride;
  uint matrix_precision;
  string frame_index_spec;
  string model_name, model_type;
  AtomicGroup model;
  MultiTrajectory mtrajA, mtrajB;
  pTraj trajectory_A, trajectory_B;
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
  SingleWorker(RealMatrix* R, vMatrix* TA, vMatrix* TB, Master* M) : _R(R), _TA(TA), _TB(TB), _M(M) { }


  SingleWorker(const SingleWorker& w) 
  {
    _R = w._R;
    _TA = w._TA;
    _TB = w._TB;
    _M = w._M;
  }
  

  void calc(const uint i) 
  {
    for (uint j=0; j<_R->cols(); ++j) 
      (*_R)(i, j) = loos::alignment::centeredRMSD((*_TA)[i], (*_TB)[j]);
  }

  void operator()() 
  {
    uint i;
    
    while (_M->workAvailable(&i))
      calc(i);
  }
  

private:
  RealMatrix* _R;
  vMatrix* _TA;
  vMatrix* _TB;
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

// just get max elt and avg like multi-rmsds.
void showStats(const RealMatrix& R) {
  uint total = R.rows() * R.cols(); //matrix is not symmetric,
  // must do calculation for each element.

  double avg = 0.0;
  double max = 0.0;
  for (uint j=0; j<R.rows(); ++j)
    for (uint i=0; i<R.cols(); ++i) {
      avg += R(j, i);
      if (R(j, i) > max)
        max = R(j, i);
    }
    
  avg /= total;
  cerr << boost::format("Max rmsd = %.4f, avg rmsd = %.4f\n") % max % avg;
}

void showFractionalStats(const RealMatrix& R, const float cutoff) {
  uint total = R.size(); 

  double avg = 0.0;
  double var = 0.0; // squared component of variance operator <x>^2 - <x^2>
  boost::tuple<uint, uint, double> max (0,0,0.0);
  //max = boost::make_tuple(0,0,0.0);
  long unsigned below_cut = 0;
  // loop over the matrix
  for (uint i=0; i<R.rows(); ++i)
    for (uint j=0; j<R.cols(); ++j) {
      if (i == j) // skip all the zeros in the diagonal
        continue;
      var += R(i, j)*R(i, j);
      avg += R(i, j);
      if (R(i, j) > boost::get<2>(max))
        boost::get<0>(max) = i;
        boost::get<1>(max) = j;
        boost::get<2>(max) = R(i, j);
      if (R(i,j) < cutoff)
        below_cut++;
    }
    
  avg /= total;
  var = var/total - avg*avg;
  cerr << boost::format("Max rmsd = %.4f between frames %d, %d, avg rmsd = %.4f, variance = %.4f, frames below %.4f = %d, total = %d\n") 
  % boost::get<2>(max) 
  % boost::get<0>(max)
  % boost::get<1>(max)
  % avg 
  % var 
  % cutoff 
  % below_cut
  % total;
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
  opts::BasicSelection* sopts = new opts::BasicSelection("name == 'backbone' && ! hydrogen");
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  verbosity = bopts->verbosity;
  report_stats = (verbosity || topts->noop);
  AtomicGroup model = topts->model;
  AtomicGroup subset = selectAtoms(model, sopts->selection);
  vector<uint> indices_A = topts->frameList(topts->trajectory_A);
  vector<uint> indices_B = topts->frameList(topts->trajectory_B);

  long mem = availableMemory();
  uint nthreads = topts->nthreads ? topts->nthreads : boost::thread::hardware_concurrency();
  
  if (verbosity > 1)
    cerr << "Using " << nthreads << " threads\n";
  
  // read in system A
  vMatrix TA = readCoords(subset, topts->trajectory_A, indices_A, verbosity > 1);
  // read in system B
  vMatrix TB = readCoords(subset, topts->trajectory_B, indices_B, verbosity > 1);

  used_memory += TA.size() * TA[0].size() * sizeof(vMatrix::value_type::value_type);   // Coords matrix for A
  used_memory += TB.size() * TB[0].size() * sizeof(vMatrix::value_type::value_type);   // Coords matrix for B
  used_memory += TA.size() * TB.size() * sizeof(RealMatrix::element_type);             // RMSDS matrix
  
  checkMemoryUsage(mem);
  centerTrajectory(TA);
  centerTrajectory(TB);

  RealMatrix M;
  if (verbosity > 1)
    cerr << "Calculating RMSD...\n";
  M = RealMatrix(TA.size(), TB.size());
  // note the 'false' here causes master to do full matrix, not just triangle.
  Master master(TA.size(), false, verbosity); 
  SingleWorker worker(&M, &TA, &TB, &master);
  Threader<SingleWorker> threads(&worker, nthreads);
  threads.join();
  if (verbosity)
    master.updateStatus();

  if (verbosity || topts->noop || topts->stats || topts->cutoff > 0){
    if (topts->cutoff > 0)
      showFractionalStats(M, topts->cutoff);
    else
      showStats(M);
  }

  if (!topts->noop) {
    cout << "# " << header << endl;
    cout << topts->trajectoryTable(topts->mtrajA);
    cout << topts->trajectoryTable(topts->mtrajB);
    cout << setprecision(topts->matrix_precision) << M;
  }

}
