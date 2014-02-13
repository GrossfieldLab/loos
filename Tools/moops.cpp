/*
  Simple molecular-order parameters for comparing CG to AA MD
  (c) 2011 Tod D. Romo, Grossfield Lab, URMC
*/


#include <loos.hpp>

using namespace std;
using namespace loos;


typedef vector<string>      vString;
typedef vector<double>      vecDbl;
typedef vector<AtomicGroup> vGroup;
typedef vector<vGroup>      vvGroup;


const double minp = 1e-3;
const double maxp = 100;
ulong nplanar = 0;
ulong ntotal = 0;


double rowStats(const RealMatrix& M, const uint row, double *std) {
  double avg = 0.0;
  for (uint i=0; i<M.cols(); ++i)
    avg += M(row, i);
  avg /= M.cols();
  
  double var = 0.0;
  for (uint i=0; i<M.cols(); ++i) {
    double d = M(row, i) - avg;
    var += d*d;
  }
  var /= (M.cols() - 1);

  *std = sqrt(var/M.cols());
  return(avg);
}


ulong calculateSize(AtomicGroup& model, const vString& names) {
  ulong n = 0;
  for (vString::const_iterator i = names.begin(); i != names.end(); ++i) {
    pTraj traj = createTrajectory(*i, model);
    n += traj->nframes();
  }
  return(n);
}


// Note: hist is really an average binned on distance
void principalComponentsOrder(dTimeSeries& order_parameters,
                              const vGroup& residues) {

  for (vGroup::const_iterator i = residues.begin(); i != residues.end(); ++i) {
    AtomicGroup residue = i->copy();
    residue.centerAtOrigin();
    residue.mergeImage();
    vector<GCoord> axes = residue.principalAxes();
    bool planar = false;

    if (axes[3].z() < minp) {
      if (nplanar == 0) {
        PDB pdb = PDB::fromAtomicGroup(residue);
        cerr << "Warning- PCA magnitudes out of bounds " << axes[3] << endl;
        cerr << pdb;
      }
      planar = true;
      ++nplanar;
    }

    double order1 = 0.5 - 1.5 * axes[1].z() * axes[1].z();
    double order2 = 0.5 - 1.5 * axes[2].z() * axes[2].z();

    order_parameters.push_back(order1);
    ++ntotal;
    if (!planar) {
      order_parameters.push_back(order2);
      ++ntotal;
    }
  }
}



vGroup extractSelections(const AtomicGroup& model, const string& selection) {
  AtomicGroup subset = selectAtoms(model, selection);
  vGroup residues = subset.splitByUniqueSegid();

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
  if (argc < 5) {
    cerr << "Usage- moops skip palm-selection model traj [traj ...] >output.asc\n";
    exit(EXIT_FAILURE);
  }

  string hdr = invocationHeader(argc, argv);
  int k = 1;


  uint skip = strtoul(argv[k++], 0, 10);
  string palm_selection = string(argv[k++]);
  AtomicGroup model = createSystem(argv[k++]);
  vGroup palms = extractSelections(model, palm_selection);

  vString traj_names;
  while (k < argc)
    traj_names.push_back(string(argv[k++]));

  cout << "# " << hdr << endl;



  // Calculate size of matrix to use...
  ulong n = calculateSize(model, traj_names);
  n -= traj_names.size() * skip;


  dTimeSeries order;

  for (vString::const_iterator i = traj_names.begin(); i != traj_names.end(); ++i) {
    dTimeSeries suborder;

    pTraj traj = createTrajectory(*i, model);
    if (skip != 0)
      traj->readFrame(skip-1);
    
    while (traj->readFrame()) {
      traj->updateGroupCoords(model);
      principalComponentsOrder(suborder, palms);
    }
    order.push_back(suborder.average());
  }

  

  cout << "Avg = " << order.average() << endl;
  cout << "Std = " << order.stdev() << endl;
  cout << boost::format("OB Data = %d out of %d (%.2f%%)\n") % nplanar % ntotal % (nplanar * 100.0 / ntotal);

  exit(EXIT_SUCCESS);
}
