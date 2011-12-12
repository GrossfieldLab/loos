/*
  contained.cpp

  Tracks the number of atoms within a blob over time...

  usage:
    contained model trajectory selection grid
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
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

#include <boost/tuple/tuple.hpp>
#include <limits>
#include <list>

#include <DensityGrid.hpp>
#include <DensityTools.hpp>


using namespace std;
using namespace loos;
using namespace loos::DensityTools;


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  opts::BasicOptions *basic_opts = new opts::BasicOptions;
  opts::BasicSelection *basic_selection = new opts::BasicSelection;
  opts::TrajectoryWithFrameIndices *basic_traj = new opts::TrajectoryWithFrameIndices;
  opts::RequiredArguments *ropts = new opts::RequiredArguments;
  ropts->addArgument("grid", "grid-name");

  opts::AggregateOptions options;
  options.add(basic_opts).add(basic_selection).add(basic_traj).add(ropts);
  if (!options.parse(argc, argv)) {
    exit(0);
  }

  AtomicGroup model = basic_traj->model;
  pTraj traj = basic_traj->trajectory;
  AtomicGroup subset = selectAtoms(model, basic_selection->selection);
  vector<uint> frames = basic_traj->frameList();

  

  cout << "# " << hdr << endl;
  cout << "# frame n\n";
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

