/*
  helix_kink.cpp

  
  Bend angle calculation over time...

  Usage:
    helix_kink selection-1 selection-2 model trajectory

  Notes:
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


string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\tAngle deviation from linear between two selections\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tGiven two selections, the angle between the first principal components\n"
    "of the two groups is calculated.  The deviation from linearity is printed\n"
    "as a function of time.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\thelix_kink 'resid <= 20' 'resid >= 25 && resid <= 44' model.psf trajectory.dcd\n"
    "Two groups are defined, the first 20 residues and residues 25 through 44.  The first\n"
    "principal component is determined for each group and the angular deviation from linear\n"
    "is printed out.\n"
    "\n";

  return(msg);
}



int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);

  if (argc != 5) {
    cerr << "Usage- helix_kink selection-1 selection-2 model trajectory\n";
    cerr << fullHelpMessage();
    exit(-1);
  }


  int k = 1;
  string pre_sel(argv[k++]);
  string post_sel(argv[k++]);
  
  AtomicGroup model = createSystem(argv[k++]);
  pTraj traj = createTrajectory(argv[k++], model);

  AtomicGroup pre = selectAtoms(model, pre_sel);
  AtomicGroup post = selectAtoms(model, post_sel);

  cout << "# " << header << endl;
  cout << boost::format("#%6s %10s     %10s %10s %10s     %10s %10s %10s\n") % "t" % "angle" % "x_0" % "y_0" % "z_0" % "x_1" % "y_1" % "z_1";

  uint t = 0;
  while (traj->readFrame()) {
    traj->updateGroupCoords(model);
    Axes preax = pre.principalAxes();
    Axes postax = post.principalAxes();

    GCoord u = preax[0];
    GCoord v = -postax[0];
    double angle = acos(u * v / (u.length() * v.length()) ) * RAD2DEG;
    if (angle > 90)
      angle = 180 - angle;
    cout << boost::format("%6d %10lf     %10lf %10lf %10lf     %10lf %10lf %10lf\n") % (t++) % angle % u[0] % u[1] % u[2] % v[0] % v[1] % v[2];
  }

}


