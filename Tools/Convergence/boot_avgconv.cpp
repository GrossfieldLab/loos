/*
  boot_avgconv.cpp
  (c) 2010 Tod D. Romo, Grossfield Lab, URMC

  Convergence of average via bootstrapping


  Usage- boot_avgconv model traj selection range repeats 0|seed
*/




#include <loos.hpp>
#include <boost/format.hpp>

using namespace loos;
using namespace std;


AtomicGroup averageSelectedSubset(const vector<AtomicGroup>& ensemble, const vector<uint>& indices) {
  AtomicGroup avg = ensemble[0].copy();
  for (AtomicGroup::iterator i = avg.begin(); i != avg.end(); ++i)
    (*i)->coords() = GCoord(0,0,0);

  uint n = avg.size();
  for (vector<uint>::const_iterator j = indices.begin(); j != indices.end(); ++j)
    for (uint i=0; i<n; ++i)
      avg[i]->coords() += ensemble[*j][i]->coords();

  for (AtomicGroup::iterator i = avg.begin(); i != avg.end(); ++i)
    (*i)->coords() /= indices.size();

  return(avg);
}




vector<uint> pickFrames(const uint nframes, const uint blocksize) {
  
  boost::uniform_int<uint> imap(0,nframes-1);
  boost::variate_generator< base_generator_type&, boost::uniform_int<uint> > rng(rng_singleton(), imap);
  vector<uint> picks;

  for (uint i=0; i<blocksize; ++i)
    picks.push_back(rng());

  return(picks);
}



int main(int argc, char *argv[]) {
  
  if (argc != 7) {
    cerr << "Usage- boot_avgconv model traj sel range nreps 0|seed\n";
    exit(0);
  }

  string hdr = invocationHeader(argc, argv);

  int k = 1;
  AtomicGroup model = createSystem(argv[k++]);
  pTraj traj = createTrajectory(argv[k++], model);
  AtomicGroup subset = selectAtoms(model, argv[k++]);
  vector<uint> sizes = parseRangeList<uint>(argv[k++]);

  uint nreps = atoi(argv[k++]);

  ulong seed = strtoul(argv[k++], 0, 10);
  if (seed == 0)
    seed = randomSeedRNG();
  else
    rng_singleton().seed(static_cast<ulong>(seed));

  cout << "# " << hdr << endl;
  cout << "# seed = " << seed << endl;
  cout << "# n\tavg\tvar\n";

  vector<AtomicGroup> ensemble;
  readTrajectory(ensemble, subset, traj);

  boost::tuple<vector<XForm>, greal, int> result = iterativeAlignment(ensemble);
  AtomicGroup avg = averageStructure(ensemble);
  
  for (vector<uint>::iterator size = sizes.begin(); size != sizes.end(); ++size) {
    TimeSeries<double> rmsds;

    for (uint n=0; n<nreps; ++n) {
      vector<uint> picks = pickFrames(ensemble.size(), *size);
      AtomicGroup sub_avg = averageSelectedSubset(ensemble, picks);
      rmsds.push_back(avg.rmsd(sub_avg));
    }

    cout << boost::format("%d\t%f\t%f\n") % *size % rmsds.average() % rmsds.variance();
  }
  
}
