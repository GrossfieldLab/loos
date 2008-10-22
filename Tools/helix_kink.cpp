/*
  bendangle.cpp

  
  (c) 2008 Tod D. Romo, Grossfield Lab
  Department of Biochemistry
  University of Rochster School of Medicine and Dentistry

  Bend angle calculation over time...

  Usage:
    bendangle selection-1 selection-2 model trajectory

  Notes:
    o Automatically appends "&& name == 'CA'" to your selection...

*/


#include <loos.hpp>

#include <boost/format.hpp>
#include <cmath>

#include "Statter.hpp"


typedef vector<GCoord> Axes;

const double RAD2DEG = 180.0 / PI;


int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);
  banal::Statter<double> stats;

  if (argc != 5) {
    cout << "Usage- bendangle selection-1 selection-2 model trajectory\n";
    exit(-1);
  }

  cout << "# " << header << endl;

  string sel1(argv[1]);
  string sel2(argv[2]);

  AtomicGroup model = loos::createSystem(argv[3]);
  pTraj ptraj = loos::createTrajectory(argv[4], model);

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
  cout << boost::format("#%10s     %10s %10s %10s     %10s %10s %10s\n") % "angle" % "x_0" % "y_0" % "z_0" % "x_1" % "y_1" % "z_1";

  while (ptraj->readFrame()) {
    ptraj->updateGroupCoords(model);
    Axes preax = pre.principalAxes();
    Axes postax = post.principalAxes();

    GCoord u = preax[0];
    GCoord v = -postax[0];
    double angle = acos(u * v / (u.length() * v.length()) ) * RAD2DEG;
    stats.push_back(angle);
    cout << boost::format("%10lf     %10lf %10lf %10lf     %10lf %10lf %10lf\n") % angle % u[0] % u[1] % u[2] % v[0] % v[1] % v[2];
  }

  cerr << stats.str() << endl;
}


