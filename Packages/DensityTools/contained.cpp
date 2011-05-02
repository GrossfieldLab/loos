/*
  contained.cpp


  (c) 2008 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry

  Tracks the number of atoms within a blob over time...

  usage:
    contained model trajectory selection grid
*/


#include <loos.hpp>

#include <boost/tuple/tuple.hpp>
#include <limits>
#include <list>

#include <DensityGrid.hpp>
#include <DensityTools.hpp>


using namespace std;
using namespace loos;
using namespace loos::DensityTools;

/// @cond INTERNAL_TOOLS

class ContainedOptions : public opts::OptionsPackage {
public:

  void addHidden(po::options_description& opts) {
    opts.add_options()
      ("grid", po::value<string>(&grid_name), "Grid name");
  }

  void addPositional(po::positional_options_description& options) {
    options.add("grid", 1);
  }

  bool check(po::variables_map& map) {
    return(!map.count("grid"));
  }

  string help() const { return("grid-name"); }
  string print() const {
    ostringstream oss;
    oss << boost::format("grid=%s") % grid_name;
    return(oss.str());
  }

  string grid_name;
};

/// @endcond


int main(int argc, char *argv[]) {
  if (argc != 5) {
    cerr << "Usage - contained model trajectory selection grid\n";
    exit(-1);
  }


  string hdr = invocationHeader(argc, argv);
  opts::BasicOptions *basic_opts = new opts::BasicOptions;
  opts::BasicSelectionOptions *basic_selection = new opts::BasicSelectionOptions;
  opts::BasicTrajectoryOptions *basic_traj = new opts::BasicTrajectoryOptions;
  ContainedOptions *my_opts = new ContainedOptions;

  opts::AggregateOptions options;
    options.add(basic_opts).add(basic_selection).add(basic_traj).add(my_opts);
  if (!options.parse(argc, argv)) {
    options.showHelp();
    exit(0);
  }

  AtomicGroup model = createSystem(basic_traj->model_name);
  pTraj traj = createTrajectory(basic_traj->traj_name, model);
  AtomicGroup subset = selectAtoms(model, basic_selection->selection);
  vector<uint> frames = opts::assignFrameIndices(traj, basic_traj->frame_index_spec, basic_traj->skip);

  

  cout << "# " << hdr << endl;
  cout << "# t n\n";
  DensityGrid<int> grid;

  ifstream ifs(my_opts->grid_name.c_str());
  if (!ifs) {
    cerr << "Error- cannot open " << my_opts->grid_name << endl;
    exit(-1);
  }
  ifs >> grid;

  for (vector<uint>::iterator i = frames.begin(); i != frames.end(); ++i) {
    traj->readFrame(*i);
    traj->updateGroupCoords(subset);

    long n = 0;
    for (AtomicGroup::iterator j = subset.begin(); j != subset.end(); ++j) {
      DensityGridpoint point = grid.gridpoint((*j)->coords());
      if (!grid.inRange(point))
	continue;
      if (grid(point) != 0)
	++n;
    }

    cout << *i << " " << n << endl;
  }
}

