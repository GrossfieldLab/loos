/*
  interdist.cpp

  (c) 2008, 2009 Tod D. Romo, Grossfield Lab
  Department of Biochemistry
  University of Rochster School of Medicine and Dentistry


  Computes distances between two selections over a trajectory...

  Assumes that the trajectory has already been aligned.

  Usage:  interdist mode pdb dcd sel1 sel2 [sel3 ...]
*/

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008,2009 Tod Romo
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
#include <limits>

using namespace std;
using namespace loos;


// Function pointer used to determine which kind of distance
// calculation we're going to do...

typedef double (*DistFPtr)(AtomicGroup&, AtomicGroup&);


// Distance between centroid of groups

double CenterDistance(AtomicGroup& u, AtomicGroup& v) {
  GCoord cu = u.centroid();
  GCoord cv = v.centroid();
  
  return(cu.distance(cv));
}


// Minimum distance between any member of group u vs any member of
// group v

double MinDistance(AtomicGroup& u, AtomicGroup& v) {
  double mind = numeric_limits<double>::max();
  AtomicGroup::iterator aj;
  for (aj = v.begin(); aj != v.end(); ++aj) {
    AtomicGroup::iterator ai;
    GCoord y = (*aj)->coords();

    for (ai = u.begin(); ai != u.end(); ++ ai) {
      double d = y.distance2((*ai)->coords());
      if (d < mind)
        mind = d;
    }
  }

  return(sqrt(mind));
}


// Maximum distance between any member of group u and any member of
// group v

double MaxDistance(AtomicGroup& u, AtomicGroup& v) {
  AtomicGroup::iterator aj;
  pAtom pv;
  
  double maxd = -1;
  for (aj = v.begin(); aj != v.end(); ++aj) {
    GCoord cv = (*aj)->coords();
    AtomicGroup::iterator ai;
    for (ai = u.begin(); ai != u.end(); ++ai) {
      double d2 = cv.distance2((*ai)->coords());
      if (d2 > maxd)
        maxd = d2;
    }
  }
  
  return (sqrt(maxd));
}




int main(int argc, char *argv[]) {
  if (argc < 6) {
    cerr << "Usage: interdist mode model trajectory sel1 sel2 [sel3 ...]\n";
    cerr << "       mode = center|min|max\n";
    exit(-1);
  }

  string header = invocationHeader(argc, argv);
  cout << "# " << header << endl;

  DistFPtr compute;
  if (strcmp(argv[1], "center") == 0)
    compute = &CenterDistance;
  else if (strcmp(argv[1], "max") == 0)
    compute = &MaxDistance;
  else if (strcmp(argv[1], "min") == 0)
    compute = &MinDistance;
  else {
    cerr << "ERROR- unknown mode '" << argv[1] << "'\n";
    cerr << "Must be either center, max, or min\n";
    exit(-1);
  }

  AtomicGroup model = createSystem(argv[2]);
  pTraj traj = createTrajectory(argv[3], model);

  AtomicGroup src = selectAtoms(model, argv[4]);

  cout << "# t ";
  
  vector<AtomicGroup> targets;
  for (int i=5; i<argc; ++i) {
    AtomicGroup trg = selectAtoms(model, argv[i]);
    targets.push_back(trg);
    cout << "d_0_" << i-4 << " ";
  }
  cout << endl;

  uint t = 0;
  while (traj->readFrame()) {

    traj->updateGroupCoords(model);

    cout << t++ << " ";

    vector<AtomicGroup>::iterator i;
    for (i = targets.begin(); i != targets.end(); ++i)
      cout << (*compute)(src, *i) << " ";
    cout << endl;
  }

}
