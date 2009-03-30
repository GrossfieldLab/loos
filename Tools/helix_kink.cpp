/*
  helix_kink.cpp

  
  (c) 2008 Tod D. Romo, Grossfield Lab
  Department of Biochemistry
  University of Rochster School of Medicine and Dentistry

  Bend angle calculation over time...

  Usage:
    helix_kink selection-1 selection-2 model trajectory

  Notes:
    o Automatically appends "&& name == 'CA'" to your selection...
    o Returns the deviation from linearity of the helix kink
      (i.e. pi - angle)
    o Angles are actually in degrees, despite the above note
      
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

#include <boost/format.hpp>
#include <cmath>

using namespace std;
using namespace loos;


typedef vector<GCoord> Axes;

const double RAD2DEG = 180.0 / PI;


int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);

  if (argc != 5) {
    cout << "Usage- helix_kink selection-1 selection-2 model trajectory\n";
    exit(-1);
  }

  cout << "# " << header << endl;

  string sel1(argv[1]);
  string sel2(argv[2]);

  AtomicGroup model = createSystem(argv[3]);
  pTraj ptraj = createTrajectory(argv[4], model);

  CAlphaSelector casel;

  Parser parsed1(sel1);
  KernelSelector ksel1(parsed1.kernel());
  AndSelector presel(ksel1, casel);
  AtomicGroup pre = model.select(presel);
  assert(pre.size() != 0 && "No atoms selected from first selection");
  cerr << boost::format("Selected %u atoms from first selection.\n") % pre.size();

  Parser parsed2(sel2);
  KernelSelector ksel2(parsed2.kernel());
  AndSelector postsel(ksel2, casel);
  AtomicGroup post = model.select(postsel);
  assert(post.size() != 0 && "No atoms selected from second selection");
  cerr << boost::format("Selected %u atoms from second selection.\n") % post.size();
  cout << boost::format("#%6s %10s     %10s %10s %10s     %10s %10s %10s\n") % "t" % "angle" % "x_0" % "y_0" % "z_0" % "x_1" % "y_1" % "z_1";

  uint t = 0;
  while (ptraj->readFrame()) {
    ptraj->updateGroupCoords(model);
    Axes preax = pre.principalAxes();
    Axes postax = post.principalAxes();

    GCoord u = preax[0];
    GCoord v = -postax[0];
    double angle = acos(u * v / (u.length() * v.length()) ) * RAD2DEG;
    cout << boost::format("%6d %10lf     %10lf %10lf %10lf     %10lf %10lf %10lf\n") % (t++) % (180.0 - angle) % u[0] % u[1] % u[2] % v[0] % v[1] % v[2];
  }

}


