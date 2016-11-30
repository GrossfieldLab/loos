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


using namespace std;
using namespace loos;


namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;


string progname;

// @cond TOOLS_INTERNAL

string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "Calculate molecular order parameters based on distance from a target\n"
    "\n"
    "DESCRIPTION\n"
    "\tDibmops is used to elucidate local effects on molecular order parameters.\n" 
    "Dibmops takes two selections, a lipopeptide (target) selection and a membrane\n"
    "lipid selection.  For each molecule in the membrane selection, the principal\n"
    "axes are determined and order parameters calculated for the 2nd and 3rd axes (as\n"
    "faux-hydrogens).  The distance to the nearest lipopeptide in the same leaflet\n"
    "is found and use to bin the order parameters.  Multiple trajectories may be\n"
    "used, in which case all trajectories are combined for binning.\n"
    "\n"
    "EXAMPLES\n"
    "\tdibmops 'resname == \"LFB\"' 'resname == \"POPC\" && name =~ \"^C2\\d+$\"' model.gro sim.xtc\n"
    "This computes a molecular order parameter for the palmitoyl chain from all POPC residues, relative\n"
    "to the LFB lipopeptides.  The default range of the histogram is [0,30) with 30 bins.\n"
    "\n"
    "\tdibmops --skip 50 --maxrad 15 --nbins 15 'resname == \"LFB\"' 'resname == \"POPC\" && name =~ \"^C2\\d+$\"' namd.psf sim1.dcd sim2.dcd\n"
    "This is the same as before, but two trajectories are used and the first 50 frames from\n"
    "each are skipped.  Additionally, the histogram range is [0,15) with 15 bins.\n"
    "\n"
    "SEE ALSO\n"
    "\tmops, order_params\n";
  

  return(msg);
}



struct ToolOptions : public opts::OptionsPackage 
{
 
  void addGeneric(po::options_description& o) 
  {
    o.add_options()
      ("skip", po::value<uint>(&skip)->default_value(0), "Skip these frames at the start of each trajectory")
      ("maxrad,R", po::value<double>(&maxrad)->default_value(30.0), "Maximum radius in membrane plane from lipopeptide")
      ("nbins,N", po::value<uint>(&nbins)->default_value(30), "Number of bins in histogram")
      ("residue", po::value<bool>(&residue_split)->default_value(false), "Force split by residue");
  
  }

  void addHidden(po::options_description& o) 
  {
    o.add_options()
      ("liposelection", po::value<string>(&liposelection), "Lipopeptide")
      ("membraneselection", po::value<string>(&membraneselection), "Membrane Lipid")
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value< vector<string> >(&traj_names), "Trajectory filenames");
  }

  void addPositional(po::positional_options_description& o) 
  {
    o.add("liposelection", 1);
    o.add("membraneselection", 1);
    o.add("model", 1);
    o.add("traj", -1);
  }
  

  bool check(po::variables_map& vm) 
  {
    return( liposelection.empty() || membraneselection.empty() || model_name.empty() || traj_names.empty() );
  }

  string help() const 
  {
    return("lipopeptide-selection membrane-lipid-selection model trajectory [trajectory ...]");
  }
  

  string print() const 
  {
    ostringstream oss;
    oss << boost::format("skip=%d, residue=%d, lipo='%s', lipid='%s', model='%s', traj='%s'")
      % skip
      % residue_split
      % liposelection
      % membraneselection
      % model_name
      % vectorAsStringWithCommas(traj_names);
    

    return(oss.str());
  }

  uint skip;
  bool residue_split;
  string membraneselection;
  string liposelection;
  string model_name;
  vector<string> traj_names;
  double maxrad;
  uint nbins;
};


typedef vector<AtomicGroup>   vecGroup;
typedef vector<double>        vecDouble;
typedef vector<ulong>         vecUlong;


enum LeafletType { UPPER, LOWER };


const double minp = 0.001;
ulong nplanar = 0;







class BinnedStatistics {
public:
  typedef vector< vector<double> >      BinType;
  typedef pair<double,double>           BinStatsType;
  
  

  BinnedStatistics(const double minval, const double maxval, const uint nbins) :
    _minval(minval),
    _nbins(nbins),
    _delta( (maxval - minval) / nbins ),
    _obdata(0),
    _npts(0),
    _bins(BinType(nbins))
  { }

  void accumulate(const double coord, const double val) {
    uint bin = (coord - _minval) / _delta;
    if (bin >= _nbins) {
      ++_obdata;
      return;
    }

    ++_npts;
    _bins[bin].push_back(val);
  }

  BinStatsType statisticsForBin(const uint bin) const {
    uint n = _bins[bin].size();
    if (n < 3)
      return(BinStatsType(0.0, 0.0));

    double mean = 0.0;
    for (uint i=0; i<n; ++i)
      mean += _bins[bin][i];
    mean /= n;

    double var = 0.0;
    for (uint i=0; i<n; ++i) {
      double d = _bins[bin][i] - mean;
      var += d*d;
    }

    return( BinStatsType(mean, sqrt(var/(n-1))/sqrt(n)) );
  }

  uint numberOfPointsForBin(const uint bin) const { return(_bins[bin].size()); }
  ulong numberOfDataPoints() const { return(_npts); }


  double binCoordinate(const uint bin) const {
    return(bin * _delta + _minval + _delta/2.0);
  }

  ulong numberOutOfBounds() const { return(_obdata); }
  uint numberOfBins() const { return(_nbins); }

private:
  double _minval;
  uint _nbins;
  double _delta;
  ulong _obdata, _npts;
  BinType _bins;
};





// @endcond 

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
void principalComponentsOrder(BinnedStatistics& phist,
                              BinnedStatistics& hist,
                              const vecGroup& molecules,
                              const vecGroup& lipopeptides) {
  vecDouble orders;

  for (vecGroup::const_iterator i = molecules.begin(); i != molecules.end(); ++i) {
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


vecGroup extractSelections(const AtomicGroup& model, const string& selection, const bool force_residues) {
  AtomicGroup subset = selectAtoms(model, selection);
  
  if (subset.empty()) {
    cerr << "Error- no atoms were selected.\n";
    exit(EXIT_FAILURE);
  }

  vecGroup residues;
  if (force_residues) {
    cerr << progname << ": Forcing split by residue\n";
    residues = subset.splitByResidue();
  } else {
    if (subset.hasBonds()) {
      cerr << progname << ": Model has connectivity.  Using this to split selection.\n";
      residues = subset.splitByMolecule();
    } else {
      residues = subset.splitByUniqueSegid();
    }
  }

  if (!force_residues && residues[0].size() == subset.size()) {
    cerr << progname << ": Either you are using a GROMACS model or you have one molecule in your selection\n";
    cerr << progname << ": If you are using GROMACS, you will want to run again with the --residue=1 option\n";
  }

  cerr << boost::format("%s: Extracted %d molecules from selection.\n") % progname % residues.size();
  return(residues);
}




int main(int argc, char *argv[]) {

  progname = string(argv[0]);
  string hdr = invocationHeader(argc, argv);
  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  ToolOptions* topts = new ToolOptions;
  
  opts::AggregateOptions options;
  options.add(bopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);
  
  uint skip = topts->skip;
  string membrane_selection = topts->membraneselection;
  string lipopeptide_selection = topts->liposelection;
  double rmax = topts->maxrad;
  uint nbins = topts->nbins;

  AtomicGroup model = createSystem(topts->model_name);
  vecGroup membrane = extractSelections(model, membrane_selection, topts->residue_split);
  vecGroup lipopeps = extractSelections(model, lipopeptide_selection, topts->residue_split);

  cerr << boost::format("Lipid selection has %d atoms per molecule and %d molecules.\n") % membrane[0].size() % membrane.size();
  cerr << boost::format("Lipopeptide selection has %d atoms per molecule and %d molecules.\n") % lipopeps[0].size() % lipopeps.size();

  BinnedStatistics lipid_phist(0.0, rmax, nbins);  // Track the dot product of the first PC with membrane
                                                   // normal (z-axis)

  BinnedStatistics lipid_hist(0.0, rmax, nbins);   // Track the fake hydrogen order parameters...


  for (vector<string>::const_iterator ci = topts->traj_names.begin(); ci != topts->traj_names.end(); ++ci) {
    pTraj traj = createTrajectory(*ci, model);
    cerr << boost::format("Processing %s ...") % *ci;
    cerr.flush();
    if (skip > 0)
      traj->readFrame(skip-1);

    while (traj->readFrame()) {
      traj->updateGroupCoords(model);


      vecGroup lip_leaf = filterByLeaflet(lipopeps, UPPER);
      if (!lip_leaf.empty()) {
          vecGroup lipid_leaf = filterByLeaflet(membrane, UPPER);
          principalComponentsOrder(lipid_phist, lipid_hist, lipid_leaf, lip_leaf);
      }

      lip_leaf = filterByLeaflet(lipopeps, LOWER);
      if (!lip_leaf.empty()) {
          vecGroup lipid_leaf = filterByLeaflet(membrane, LOWER);
          principalComponentsOrder(lipid_phist, lipid_hist, lipid_leaf, lip_leaf);
      }
    }

    cerr << " done\n";
  }



  cerr << boost::format("Lipid histogram had %d points with %d out-of-bounds\n") % lipid_hist.numberOfDataPoints() % lipid_hist.numberOutOfBounds();

  cout << "# " << hdr << endl;
  cout << "# Lipid total = " << lipid_hist.numberOfDataPoints() << endl;
  cout << "# d\tLipid-n\tLipid-avg\tLipid-stderr\tLipid-1stPC\n";


  for (uint i=0; i<nbins; ++i) {
    BinnedStatistics::BinStatsType stats = lipid_hist.statisticsForBin(i);
    cout << lipid_hist.binCoordinate(i) << "\t" << lipid_hist.numberOfPointsForBin(i) << "\t" << stats.first << "\t" << stats.second << "\t";

    stats = lipid_phist.statisticsForBin(i);
    cout << stats.first << "\n";
  }

  if (nplanar != 0)
    cerr << boost::format("Warning- there were %d planar lipids found\n") % nplanar;

}
