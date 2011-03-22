/*
  water-inside.cpp


  (c) 2008 Tod D. Romo,
      Grossfield Lab,
      University of Rochester Medical and Dental School


  Applies a given set of criteria to determine whether or not a water
  is inside a protein.  A matrix is then built up where each column
  represents a timepoint in the trajectory and each row is the
  internal water state (i.e. 1 = water is inside, 0 = water is not
  inside)

  Also tracks the volume of the probe region (i.e. what's defined as
  inside, if possible) and writes out a list of atomids that describe
  which atoms go with which rows of the matrix.

*/

#include <loos.hpp>



#include <loos.hpp>
#include <cmath>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <DensityGrid.hpp>
#include <internal-water-filter.hpp>
#include <OptionsFramework.hpp>
#include <DensityOptions.hpp>

using namespace std;
using namespace loos;
using namespace loos::DensityTools;
namespace opts = loos::DensityTools::OptionsFramework;



typedef Math::Matrix<int, Math::ColMajor>    Matrix;


// Globals...  BOOO!
WaterFilterBase* filter_func;
string water_string, prot_string, model_name, traj_name, prefix;
DensityGrid<int> the_grid;



    

void writeAtomIds(const string& fname, const AtomicGroup& grp, const string& hdr) {
  ofstream ofs(fname.c_str());
  AtomicGroup::const_iterator ci;

  uint t = 0;
  ofs << "# " << hdr << endl;
  ofs << "# i\tatomid(i)\tresidue(i)\n";
  for (ci = grp.begin(); ci != grp.end(); ++ci) {
    ofs << boost::format("%u\t%d\t%s-%d\n") % t++ % (*ci)->id() % (*ci)->name() % (*ci)->resid();
  }
}


int main(int argc, char *argv[]) {
  string hdr = loos::invocationHeader(argc, argv);

  opts::BasicOptions* basopts = new opts::BasicOptions;
  opts::OutputPrefixOptions* prefopts = new opts::OutputPrefixOptions;
  opts::BasicTrajectoryOptions* trajopts = new opts::BasicTrajectoryOptions;
  opts::BasicWaterOptions* watopts = new opts::BasicWaterOptions;

  opts::AggregateOptions options;
  options.addOptions(basopts).addOptions(prefopts).addOptions(trajopts).addOptions(watopts);
  if (!options.parseOptions(argc, argv))
    exit(-1);

  // Copy globals back
  model_name = trajopts->model_name;
  traj_name = trajopts->traj_name;
  prot_string = watopts->prot_string;
  water_string = watopts->water_string;
  filter_func = watopts->filter_func;
  the_grid = watopts->the_grid;
  prefix = prefopts->prefix;


  AtomicGroup model = loos::createSystem(model_name);
  pTraj traj = loos::createTrajectory(traj_name, model);
  vector<uint> frames = assignFrameIndices(traj, trajopts->frame_index_spec, trajopts->skip);

  AtomicGroup subset = loos::selectAtoms(model, prot_string);
  AtomicGroup waters = loos::selectAtoms(model, water_string);

  uint m = waters.size();
  uint n = traj->nframes();
  Matrix M(m,n);
  Math::Matrix<double> V(n, 1);
  cerr << boost::format("Water matrix is %d x %d.\n") % m % n;

  uint i = 0;
  cerr << "Processing - ";

  for (vector<uint>::iterator t = frames.begin(); t != frames.end(); ++t) {
    if (i % 100 == 0)
      cerr << ".";

    traj->readFrame(*t);
    traj->updateGroupCoords(model);

    vector<int> mask = filter_func->filter(waters, subset);
    if (mask.size() != m) {
      cerr << boost::format("ERROR - returned mask has size %u but expected %u.\n") % mask.size() % m;
      exit(-10);
    }

    for (uint j=0; j<m; ++j)
      M(j,i) = mask[j];

    V(i,0) = filter_func->volume();
    ++i;
  }

  cerr << " done\n";
  writeAsciiMatrix(prefix + ".asc", M, hdr);
  writeAsciiMatrix(prefix + ".vol", V, hdr);
  writeAtomIds(prefix + ".atoms", waters, hdr);
}
