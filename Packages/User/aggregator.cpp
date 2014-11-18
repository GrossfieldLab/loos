// (c) 2012 Tod D. Romo, Grossfield Lab, URMC

#include <loos.hpp>

#include <list>

using namespace std;
using namespace loos;


typedef vector<AtomicGroup>    vGroup;
typedef list<AtomicGroup>      lGroup;
typedef vector<uint>           vUint;


bool inContact(const AtomicGroup& mol1, const AtomicGroup& mol2, const double cutoff, const uint ncontacts) {
  uint count = 0;
  double cutoff2 = cutoff * cutoff;
  GCoord box = mol1.periodicBox();

  for (AtomicGroup::const_iterator j = mol1.begin(); j != mol1.end(); ++j) {
    GCoord c = (*j)->coords();
    for (AtomicGroup::const_iterator i = mol2.begin(); i != mol2.end(); ++i) {
      double d = c.distance2( (*i)->coords(), box );
      if (d <= cutoff2)
        if (++count >= ncontacts)
          return(true);
    }
  }

  return(false);
}



double avgClusterSize(const lGroup& clusters) {
  double avg = 0.0;

  for (lGroup::const_iterator i = clusters.begin(); i != clusters.end(); ++i)
    avg += i->size();

  avg /= clusters.size();
  return(avg);
}


double avgRadius(const lGroup& clusters) {
  double avg = 0.0;

  for (lGroup::const_iterator i = clusters.begin(); i != clusters.end(); ++i)
    avg += i->radius();

  avg /= clusters.size();
  return(avg);
}



int main(int argc, char *argv[]) {

  if (argc != 6) {
    cerr << "Usage- aggregator model traj selection #-of-contacts contact-distance\n";
    exit(-1);
  }
  
  string hdr = invocationHeader(argc, argv);
  int k = 1;
  AtomicGroup model = createSystem(argv[k++]);
  pTraj traj = createTrajectory(argv[k++], model);
  if (!traj->hasPeriodicBox()) {
    cerr << "Error- trajectory has no periodic boundary information.\n";
    exit(-2);
  }

  string selection(argv[k++]);
  uint ncontacts = strtoul(argv[k++], 0, 10);
  double dcutoff = strtod(argv[k++], 0);

  AtomicGroup subset = selectAtoms(model, selection);
  vGroup molecules = subset.splitByMolecule();
  if (molecules.size() <= 1) {
    cerr << "Error- you need at least two molecules.\n";
    exit(-2);
  }

  cout << "# " << hdr << endl;
  cout << boost::format("# Found %d molecules\n") % molecules.size();
  cout << "# t number-of-clusters\tavg-atoms-per-cluster\tavg-radius-per-cluster\n";

  uint t = 0;
  while (traj->readFrame()) {
    traj->updateGroupCoords(model);

    lGroup clusters;
    clusters.push_back(molecules[0]);

    for (vGroup::iterator mol = molecules.begin()+1; mol != molecules.end(); ++mol) {
      vector<lGroup::iterator> contacts;

      for (lGroup::iterator i = clusters.begin(); i != clusters.end(); ++i)
        if (inContact(*mol, *i, dcutoff, ncontacts))
          contacts.push_back(i);

      if (contacts.empty())
        clusters.push_back(*mol);
      else {
        AtomicGroup new_cluster = *mol;
        for (vector<lGroup::iterator>::iterator i = contacts.begin(); i != contacts.end(); ++i) {
          new_cluster.append(**i);
          clusters.erase(*i);
        }
        clusters.push_back(new_cluster);
      }

    }

    cout << boost::format("%d\t%d\t%f\t%f\n")
      % (t++)
      % clusters.size()
      % avgClusterSize(clusters)
      % avgRadius(clusters);
  }
}
