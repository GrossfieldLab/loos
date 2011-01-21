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
  
  if (argc < 5 || argc > 6) {
    cerr << "Usage- block_avgconv model traj sel range [1 = do not align trajectory]\n";
    exit(0);
  }

  string hdr = invocationHeader(argc, argv);

  int k = 1;
  AtomicGroup model = createSystem(argv[k++]);
  pTraj traj = createTrajectory(argv[k++], model);
  AtomicGroup subset = selectAtoms(model, argv[k++]);
  vector<uint> sizes = parseRangeList<uint>(argv[k++]);
  bool do_align = argc != 6;

  cout << "# " << hdr << endl;
  cout << "# n\tavg\tvar\tblocks\tstderr\n";

  vector<AtomicGroup> ensemble;
  cerr << "Reading trajectory...\n";
  readTrajectory(ensemble, subset, traj);

  if (do_align) {
    cerr << "Aligning trajectory...\n";
    boost::tuple<vector<XForm>, greal, int> result = iterativeAlignment(ensemble);
  } else
    cerr << "Trajectory is already aligned!\n";

  cerr << "Processing- ";
  for (uint block = 0; block < sizes.size(); ++block) {
    if (block % 50)
      cerr << ".";

    uint blocksize = sizes[block];

    vector<AtomicGroup> averages;
    for (uint i=0; i<ensemble.size() - blocksize; i += blocksize) {

      vector<uint> indices(blocksize);
      for (uint j=0; j<blocksize; ++j)
        indices[j] = i+j;
      
      averages.push_back(averageSelectedSubset(ensemble, indices));
    }
    
    TimeSeries<double> rmsds;
    for (uint j=0; j<averages.size() - 1; ++j)
      for (uint i=j+1; i<averages.size(); ++i) {
        AtomicGroup left = averages[j];
        AtomicGroup right = averages[i];
        left.alignOnto(right);
        rmsds.push_back(left.rmsd(right));
      }
        

    double v = rmsds.variance();
    uint n = averages.size();
    cout << boost::format("%d\t%f\t%f\t%d\t%f\n") % blocksize % rmsds.average() % v % n % sqrt(v/n);
  }
  
  cerr << "\nDone!\n";
}
