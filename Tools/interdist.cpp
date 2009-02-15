/*
  interdist.cpp

  (c) 2008, 2009 Tod D. Romo, Grossfield Lab
  Department of Biochemistry
  University of Rochster School of Medicine and Dentistry


  Computes distances between two selections over a trajectory...

  Assumes that the trajectory has already been aligned.

  Usage:  interdist mode pdb dcd sel1 sel2 [sel3 ...]
*/




#include <loos.hpp>
#include <getopt.h>

using namespace std;
using namespace loos;

typedef double (*DistFPtr)(AtomicGroup&, AtomicGroup&);

double CenterDistance(AtomicGroup& u, AtomicGroup& v) {
  GCoord cu = u.centroid();
  GCoord cv = v.centroid();
  
  return(cu.distance(cv));
}


double MinDistance(AtomicGroup& u, AtomicGroup& v) {
  AtomicGroup::Iterator viter(v);
  pAtom pv;
    
  double mind = 1e30;
  while (pv = viter()) {
    GCoord cv = pv->coords();
    AtomicGroup::Iterator uiter(u);
    pAtom pu;
    while (pu = uiter()) {
      GCoord cu = pu->coords();
      double d2 = cv.distance2(cu);
      if (d2 < mind)
        mind = d2;
    }
  }
  
  return (sqrt(mind));
}


double MaxDistance(AtomicGroup& u, AtomicGroup& v) {
  AtomicGroup::Iterator viter(v);
  pAtom pv;
  
  double maxd = -1;
  while (pv = viter()) {
    GCoord cv = pv->coords();
    AtomicGroup::Iterator uiter(u);
    pAtom pu;
    while (pu = uiter()) {
      GCoord cu = pu->coords();
      double d2 = cv.distance2(cu);
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
