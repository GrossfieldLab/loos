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
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/tuple/tuple.hpp>
#include <limits>
#include <list>

#include <DensityGrid.hpp>


using namespace std;
using namespace loos;
using namespace loos::DensityTools;




int main(int argc, char *argv[]) {
  if (argc != 5) {
    cerr << "Usage - contained model trajectory selection grid\n";
    exit(-1);
  }


  string hdr = invocationHeader(argc, argv);
  cout << "# " << hdr << endl;
  cout << "# t n\n";

  AtomicGroup model = createSystem(argv[1]);
  pTraj traj = createTrajectory(argv[2], model);

  AtomicGroup subset = selectAtoms(model, argv[3]);

  ifstream ifs(argv[4]);
  if (!ifs) {
    cerr << "Error - cannot open " << argv[4] << " for reading.\n";
    exit(-1);
  }
  DensityGrid<int> grid;

  ifs >> grid;

  uint t = 0;

  while (traj->readFrame()) {
    traj->updateGroupCoords(subset);

    long n = 0;

    AtomicGroup::Iterator ai(subset);
    pAtom atom;

    while (atom = ai()) {
      DensityGridpoint point = grid.gridpoint(atom->coords());
      if (!grid.inRange(point))
	continue;
      if (grid(point) != 0)
	++n;
    }

    cout << t++ << " " << n << endl;
  }
}
