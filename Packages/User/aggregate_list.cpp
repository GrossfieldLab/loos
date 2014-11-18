/*
 *
 * Hacked up Tod's aggregator tool - returns lists of clusters and
 * what molecules make them up. Useful for manual curation of
 * clusters...
 *
 */

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

  if (argc != 7) {
    cerr << "Usage- aggregator model traj selection #-of-contacts contact-distance endframe \n";
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

  uint end = strtoul(argv[k++],0, 10);

  AtomicGroup subset = selectAtoms(model, selection);
  vGroup molecules = subset.splitByMolecule();
  if (molecules.size() <= 1) {
    cerr << "Error- you need at least two molecules.\n";
    exit(-2);
  }

  cout << "# " << hdr << endl;
  cout << boost::format("# Found %d molecules\n") % molecules.size();
  cout << "# t number-of-clusters\tavg-atoms-per-cluster\tavg-radius-per-cluster\n";

  uint frame = 0;
  uint t = 0;
  while (traj->readFrame() && frame<end) {
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

    // Loop over the clusters and perform the important stuff!!!
    
    for (lGroup::iterator clust = clusters.begin(); clust != clusters.end(); ++clust) {

      vGroup clustMol = clust->splitByMolecule();
      string ends = "";
      string cluster = "";


      for (vGroup::iterator m = clustMol.begin(); m != clustMol.end(); ++m) {

        int contactCount = 0;


        for (vGroup::iterator n = clustMol.begin(); n != clustMol.end(); ++n) {

            if ( m != n ){


                if (inContact(*m, *n, dcutoff, ncontacts))
                    contactCount++;

            }

        }

        if (contactCount == 1 && clustMol.size() > 1) {
            ends = ends + m[0].getAtom(0)->segid().at(2) + m[0].getAtom(0)->segid().at(3) + " ";
        }
        

        cluster = cluster + m[0].getAtom(0)->segid().at(2) + m[0].getAtom(0)->segid().at(3) + "|" ;

      }

      if (ends.empty()){

        ends = "NONE NONE";

      }

      cout << frame << "\tSize:" << clustMol.size() << "\t";
      cout << ends << "\t" << cluster.substr(0, cluster.size()-1) << "\t";
      cout << endl;

    }

    //cout << boost::format("%d\t%d\t%f\t%f\n")
    //  % (t++)
    //  % clusters.size()
    //  % avgClusterSize(clusters)
    //  % avgRadius(clusters);
  
  frame++;
  }
}
