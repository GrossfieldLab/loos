/*
  DIstance Based Molecular Order Parameters
*/

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2011, Tod D. Romo and Alan Grossfield
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
#include <shist.hpp>


using namespace std;
using namespace loos;

typedef vector<AtomicGroup>   vecGroup;
typedef vector<double>        vecDouble;
typedef vector<ulong>         vecUlong;


enum LeafletType { UPPER, LOWER };


const double minp = 0.001;
const double maxp = 100;
ulong nplanar = 0;

// Only look at x,y-plane distance
double centroidDistance(const AtomicGroup& a, const AtomicGroup& b) {
  GCoord ac = a.centroid();
  ac.z() = 0.0;

  GCoord bc = b.centroid();
  bc.z() = 0.0;
  return(ac.distance(bc));
}

// Uses x,y-plane distance between centroids
double minDistanceToSet(const AtomicGroup& a, const vecGroup& set) {

  double mind = numeric_limits<double>::max();
  for (vecGroup::const_iterator i = set.begin(); i != set.end(); ++i) {
    double d = centroidDistance(a, *i);
    if (d < mind)
      mind = d;
  }

  return(mind);
}




// Note: hist is really an average binned on distance
void principalComponentsOrder(BinnedStatistics<double>& phist,
                              BinnedStatistics<double>& hist,
                              const vecGroup& residues,
                              const vecGroup& lipopeptides) {
  vecDouble orders;

  for (vecGroup::const_iterator i = residues.begin(); i != residues.end(); ++i) {
    vector<GCoord> axes = i->principalAxes();
    bool planar = false;

    if (axes[3].z() < minp) {
      if (nplanar == 0) {
        PDB pdb = PDB::fromAtomicGroup(*i);
        cerr << "Warning- PCA magnitudes out of bounds " << axes[3] << endl;
        cerr << pdb;
      }
      planar = true;
      ++nplanar;
    }


    double order1 = 0.5 - 1.5 * axes[1].z() * axes[1].z();
    double order2 = 0.5 - 1.5 * axes[2].z() * axes[2].z();

    GCoord lipid_center = i->centroid();
    lipid_center.z() = 0.0;

    for (vecGroup::const_iterator j = lipopeptides.begin(); j != lipopeptides.end(); ++j) {
      GCoord c = j->centroid();
      c.z() = 0.0;
      double d = lipid_center.distance(c);
      phist.accumulate(d, abs(axes[0].z()));
      
      hist.accumulate(d, order1);
      if (!planar)
        hist.accumulate(d, order2);
    }
    
  }
}


vecGroup filterByLeaflet(const vecGroup& ensemble, const LeafletType leaflet = UPPER) {
  vecGroup side;

  for (vecGroup::const_iterator i = ensemble.begin(); i != ensemble.end(); ++i) {
    GCoord c = i->centroid();
    if (leaflet == UPPER && c.z() > 0)
      side.push_back(*i);
    else if (leaflet == LOWER && c.z() < 0)
      side.push_back(*i);
  }

  return(side);
}


vecGroup extractSelections(const AtomicGroup& model, const string& selection) {
  AtomicGroup subset = selectAtoms(model, selection);
  vecGroup residues = subset.splitByUniqueSegid();

  if (residues.empty()) {
    cerr << boost::format("ERROR- could not split group using selection '%s'\n") % selection;
    exit(EXIT_FAILURE);
  }
  
  // Autodetect whether we should use segid or residue to split...
  if (residues[0].size() == subset.size()) {
    cerr << "WARNING- apparent GROMACS source data...switching to splitByResidue() mode\n";
    residues = subset.splitByResidue();
  }
  return(residues);
}




int main(int argc, char *argv[]) {

  if (argc < 9) {
    cerr << "Usage- dibaprop skip palm-selection oleo-selection lipopeptide-selection xmax xbins  model traj [traj...]\n";
    exit(EXIT_FAILURE);
  }


  string hdr = invocationHeader(argc, argv);
  int k = 1;

  uint skip = strtoul(argv[k++], 0, 10);
  string palm_selection = string(argv[k++]);
  string oleo_selection = string(argv[k++]);
  string lipopeptide_selection = string(argv[k++]);
  double xmax = strtod(argv[k++], 0);
  uint xbins = strtoul(argv[k++], 0, 10);

  AtomicGroup model = createSystem(argv[k++]);
  vecGroup palms = extractSelections(model, palm_selection);
  vecGroup oleos = extractSelections(model, oleo_selection);
  vecGroup lipopeps = extractSelections(model, lipopeptide_selection);

  cerr << boost::format("Palm selection has %d atoms per residue and %d residues.\n") % palms[0].size() % palms.size();
  cerr << boost::format("Oleo selection has %d atoms per residue and %d residues.\n") % oleos[0].size() % oleos.size();
  cerr << boost::format("Lipopeptide selection has %d atoms per residue and %d residues.\n") % lipopeps[0].size() % lipopeps.size();

  BinnedStatistics<double> palm_phist(0.0, xmax, xbins);
  BinnedStatistics<double> palm_hist(0.0, xmax, xbins);

  BinnedStatistics<double> oleo_phist(0.0, xmax, xbins);
  BinnedStatistics<double> oleo_hist(0.0, xmax, xbins);

  while (k < argc) {
    pTraj traj = createTrajectory(argv[k++], model);
    cerr << boost::format("Processing %s ...") % argv[k-1];
    cerr.flush();
    if (skip > 0)
      traj->readFrame(skip-1);

    while (traj->readFrame()) {
      traj->updateGroupCoords(model);


      vecGroup lip_leaf = filterByLeaflet(lipopeps, UPPER);
      if (!lip_leaf.empty()) {
          vecGroup palm_leaf = filterByLeaflet(palms, UPPER);
          vecGroup oleo_leaf = filterByLeaflet(oleos, UPPER);
          principalComponentsOrder(palm_phist, palm_hist, palm_leaf, lip_leaf);
          principalComponentsOrder(oleo_phist, oleo_hist, oleo_leaf, lip_leaf);
      }

      lip_leaf = filterByLeaflet(lipopeps, LOWER);
      if (!lip_leaf.empty()) {
          vecGroup palm_leaf = filterByLeaflet(palms, LOWER);
          vecGroup oleo_leaf = filterByLeaflet(oleos, LOWER);
          principalComponentsOrder(palm_phist, palm_hist, palm_leaf, lip_leaf);
          principalComponentsOrder(oleo_phist, oleo_hist, oleo_leaf, lip_leaf);
      }
    }

    cerr << " done\n";
  }



  cerr << boost::format("Palm histogram had %d points with %d out-of-bounds\n") % palm_hist.numberOfDataPoints() % palm_hist.numberOutOfBounds();
  cerr << boost::format("Oleo histogram had %d points with %d out-of-bounds\n") % oleo_hist.numberOfDataPoints() % oleo_hist.numberOutOfBounds();

  cout << "# " << hdr << endl;
  cout << "# Palm total = " << palm_hist.numberOfDataPoints() << endl;
  cout << "# Oleo total = " << oleo_hist.numberOfDataPoints() << endl;
  cout << "# d\tPalm-n\tPalm-avg\tPalm-stderr\tOleo-n\tOleo-avg\tOleo-stderr\tPalm-1stPC\tOloe-1stPC\n";


  for (uint i=0; i<xbins; ++i) {
    BinnedStatistics<double>::BinStatsType stats = palm_hist.statisticsForBin(i);
    cout << palm_hist.binCoordinate(i) << "\t" << palm_hist.numberOfPointsForBin(i) << "\t" << stats.first << "\t" << stats.second << "\t";

    stats = oleo_hist.statisticsForBin(i);
    cout << "\t" << oleo_hist.numberOfPointsForBin(i) << "\t" << stats.first << "\t" << stats.second << "\t";

    stats = palm_phist.statisticsForBin(i);
    cout << stats.first << "\t";
    stats = oleo_phist.statisticsForBin(i);
    cout << stats.first << endl;
  }

  if (nplanar != 0)
    cerr << boost::format("Warning- there were %d planar lipids found\n") % nplanar;

}
