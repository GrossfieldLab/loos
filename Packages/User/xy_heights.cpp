/*
  Position along the z-axis for a selection as a function of distance
  from another selection.

  I used it to calculated bilayer curvature induced by lipopeptides
  (essentially phosphate heights for each leaflet as a function of distance
  from a lipopeptide)

  Essentially hacked up from Tod's dibaprop tool...
  (c) 2011 Tod D. Romo, Grossfield Lab, URMC
*/



#include <loos.hpp>
#include <shist.hpp>


using namespace std;
using namespace loos;

typedef vector<AtomicGroup>   vecGroup;
typedef vector<double>        vecDouble;
typedef vector<ulong>         vecUlong;


enum LeafletType { UPPER, LOWER };


const double minp = 0.0001;
const double maxp = 1000;

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
void heightMap(BinnedStatistics<double>& hist,
                              const vecGroup& residues,
                              const vecGroup& lipopeptides) {

  for (vecGroup::const_iterator i = residues.begin(); i != residues.end(); ++i) {

    double height = (i->centroid()).z();

    GCoord lipid_center = i->centroid();
    lipid_center.z() = 0.0;

    for (vecGroup::const_iterator j = lipopeptides.begin(); j != lipopeptides.end(); ++j) {
      GCoord c = j->centroid();
      c.z() = 0.0;
      double d = lipid_center.distance(c);
      
      hist.accumulate(d, height);
    
    
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
  vecGroup residues = subset.splitByMolecule();

  if (residues.empty()) {
    cerr << boost::format("ERROR- could not split group using selection '%s'\n") % selection;
    exit(EXIT_FAILURE);
  }
  
  // Autodetect whether we should use segid or residue to split...
  if (residues[0].size() == subset.size()) {
    cerr << "WARNING- apparent GROMACS source data...switching to splitByResidue() mode\n";
    residues = subset.splitByMolecule();
  }
  return(residues);
}




int main(int argc, char *argv[]) {

  if (argc < 8) {
    cerr << "Usage- xy_heights skip lipid-selection lipopeptide-selection xmax xbins model traj [traj...]\n";
    exit(EXIT_FAILURE);
  }


  string hdr = invocationHeader(argc, argv);
  int k = 1;

  uint skip = strtoul(argv[k++], 0, 10);
  string lipid_selection = string(argv[k++]);
  string lipopeptide_selection = string(argv[k++]);
  double xmax = strtod(argv[k++], 0);
  uint xbins = strtoul(argv[k++], 0, 10);

  AtomicGroup model = createSystem(argv[k++]);
  vecGroup lipids = extractSelections(model, lipid_selection);
  vecGroup lipopeps = extractSelections(model, lipopeptide_selection);

  cerr << boost::format("Lipid selection has %d atoms per residue and %d residues.\n") % lipids[0].size() % lipids.size();
  cerr << boost::format("Lipopeptide selection has %d atoms per residue and %d residues.\n") % lipopeps[0].size() % lipopeps.size();

  BinnedStatistics<double> lipid_upper_hist(0.0, xmax, xbins);
  BinnedStatistics<double> lipid_lower_hist(0.0, xmax, xbins);

  while (k < argc) {
    pTraj traj = createTrajectory(argv[k++], model);
    cerr << boost::format("Processing %s ...") % argv[k-1];
    cerr.flush();
    if (skip > 0)
      traj->readFrame(skip-1);

    while (traj->readFrame()) {
      traj->updateGroupCoords(model);


          vecGroup lipid_leaf = filterByLeaflet(lipids, UPPER);
          heightMap(lipid_upper_hist, lipid_leaf, lipopeps);

          lipid_leaf = filterByLeaflet(lipids, LOWER);
          heightMap(lipid_lower_hist, lipid_leaf, lipopeps);
      }
    }

    cerr << " done\n";


  cout << "# " << hdr << endl;
  cout << "# Upper lipid total = " << lipid_upper_hist.numberOfDataPoints() << endl;
  cout << "# Lower lipid total = " << lipid_lower_hist.numberOfDataPoints() << endl;
  cout << "# d\tUpper\tavg\tstderror\tLower\tavg\tstderror\n";


  for (uint i=0; i<xbins; ++i) {
    BinnedStatistics<double>::BinStatsType stats = lipid_upper_hist.statisticsForBin(i);
    BinnedStatistics<double>::BinStatsType stats_low = lipid_lower_hist.statisticsForBin(i);

    cout << lipid_upper_hist.binCoordinate(i) << "\t" << lipid_upper_hist.numberOfPointsForBin(i) << "\t" << stats.first << "\t" << stats.second << "\t";
    cout << lipid_lower_hist.binCoordinate(i) << "\t" << lipid_lower_hist.numberOfPointsForBin(i) << "\t" << stats_low.first << "\t" << stats_low.second << "\t";
    cout << endl;
  }


}
