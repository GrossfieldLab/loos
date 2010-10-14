/*
  avgconv_drift
  (c) 2010 Tod D. Romo, Grossfield Lab, URMC

  Usage- avgconv model traj selection range
*/

#include <loos.hpp>
#include <boost/format.hpp>

using namespace std;
using namespace loos;


bool locally_optimal = false;


AtomicGroup calcAverage(const vector<AtomicGroup>& ensemble, const uint size) {

  vector<AtomicGroup> subsample(size);
  for (uint i=0; i<size; ++i)
    subsample[i] = ensemble[i];

  if (locally_optimal)
    (void)iterativeAlignment(subsample);

  AtomicGroup avg = averageStructure(subsample);
  return(avg);
}



int main(int argc, char *argv[]) {
  if (argc < 5 || argc > 6) {
    cout << "Usage- avgconv model traj selection step [1 = local optimal avg]\n";
    exit(-1);
  }

  string hdr = invocationHeader(argc, argv);
  cout << "# " << hdr << endl;
  cout << "# n\trmsd\n";

  int k = 1;
  AtomicGroup model = createSystem(argv[k++]);
  pTraj traj = createTrajectory(argv[k++], model);
  string sel = string(argv[k++]);
  uint step = atoi(argv[k++]);
  locally_optimal =  (argv[k++][0] == '1');
  
  AtomicGroup subset = selectAtoms(model, sel);
  cout << boost::format("# Subset has %d atoms\n") % subset.size();
  vector<AtomicGroup> ensemble;
  readTrajectory(ensemble, subset, traj);
  cout << boost::format("# Trajectory has %d frames\n") % ensemble.size();

  if (!locally_optimal) {
    boost::tuple<vector<XForm>, greal, int> result = iterativeAlignment(ensemble);
    cout << boost::format("# Iterative alignment converged to RMSD of %g with %d iterations\n") % boost::get<1>(result) % boost::get<2>(result);
  }

  AtomicGroup preceding = calcAverage(ensemble, step);
  for (uint n=2*step; n<ensemble.size(); n += step) {
    AtomicGroup avg = calcAverage(ensemble, n);
    double rmsd = preceding.rmsd(avg);

    cout << boost::format("%d\t%f\n") % n % rmsd;

    preceding = avg;
  }

}

