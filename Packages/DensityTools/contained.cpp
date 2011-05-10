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


int main(int argc, char *argv[]) {
  if (argc != 5) {
    cerr << "Usage - contained model trajectory selection grid\n";
    exit(-1);
  }


  string hdr = invocationHeader(argc, argv);
  opts::BasicOptions *basic_opts = new opts::BasicOptions;
  opts::BasicSelectionOptions *basic_selection = new opts::BasicSelectionOptions;
  opts::BasicTrajectoryOptions *basic_traj = new opts::BasicTrajectoryOptions;
  opts::RequiredOptions *ropts = new opts::RequiredOptions;
  ropts->addOption("grid", "grid-name");

  opts::AggregateOptions options;
  options.add(basic_opts).add(basic_selection).add(basic_traj).add(ropts);
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

  string grid_name = ropts->value("grid");
  ifstream ifs(grid_name.c_str());
  if (!ifs) {
    cerr << "Error- cannot open " << grid_name << endl;
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

