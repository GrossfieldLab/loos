/*
  block_avgconv.cpp
  (c) 2010 Tod D. Romo, Grossfield Lab, URMC

  Convergence of average via block averaging

  Usage- block_avgconv model traj selection range
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




int main(int argc, char *argv[]) {
  
  if (argc != 5) {
    cerr << "Usage- block_avgconv model traj sel range\n";
    exit(0);
  }

  string hdr = invocationHeader(argc, argv);

  int k = 1;
  AtomicGroup model = createSystem(argv[k++]);
  pTraj traj = createTrajectory(argv[k++], model);
  AtomicGroup subset = selectAtoms(model, argv[k++]);
  vector<uint> sizes = parseRangeList<uint>(argv[k++]);

  cout << "# " << hdr << endl;
  cout << "# n\tavg\tvar\n";

  vector<AtomicGroup> ensemble;
  readTrajectory(ensemble, subset, traj);

  boost::tuple<vector<XForm>, greal, int> result = iterativeAlignment(ensemble);
  AtomicGroup avg = averageStructure(ensemble);
  
  for (vector<uint>::iterator size = sizes.begin(); size != sizes.end(); ++size) {
    TimeSeries<double> rmsds;

    for (uint i=0; i<ensemble.size() - *size; i += *size) {
      vector<uint> indices(*size);
      for (uint j=0; j<*size; ++j)
        indices[j] = i+j;
      
      AtomicGroup sub_avg = averageSelectedSubset(ensemble, indices);
      rmsds.push_back(avg.rmsd(sub_avg));
    }

    cout << boost::format("%d\t%f\t%f\n") % *size % rmsds.average() % rmsds.variance();
  }
  
}
