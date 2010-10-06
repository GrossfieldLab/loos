/*
  
  (c) 2009,2010 Tod D. Romo, Grossfield Lab, URMC...

  Perform a block-overlap in comparison to a full PCA


  usage- bcom model trajectory selection modes block-list

  Notes:
    o if modes = 0, then use all available modes
*/


#include <loos.hpp>
#include "bcomlib.hpp"

using namespace std;
using namespace loos;
using namespace Convergence;



typedef vector<AtomicGroup>                               vGroup;
typedef boost::tuple<RealMatrix, RealMatrix, RealMatrix>  SVDResult;


struct Datum {
  Datum(const double avg, const double var, const double power) : avg_coverlap(avg),
                                                                  var_coverlap(var),
                                                                  avg_power(power) { }


  double avg_coverlap;
  double var_coverlap;
  double avg_power;
};



const bool length_normalize = true;




vGroup subgroup(const vGroup& A, const uint a, const uint b) {
  vGroup B;

  for (uint i=a; i<b; ++i)
    B.push_back(A[i]);

  return(B);
}


double sum(const RealMatrix& v) {
  double s = 0.0;
  for (uint j=0; j<v.rows(); ++j)
    s += v[j];

  return(s);
}


template<class ExtractPolicy>
Datum blocker(const RealMatrix& Ua, const RealMatrix sa, vGroup& ensemble, const uint blocksize, ExtractPolicy& policy) {


  vector<double> coverlaps;
  vector<double> powers;

  double sa_sum = sum(sa);

  for (uint i=0; i<ensemble.size() - blocksize; i += blocksize) {
    vGroup subset = subgroup(ensemble, i, i+blocksize);
    boost::tuple<RealMatrix, RealMatrix> pca_result = pca(subset, policy);
    RealMatrix s = boost::get<0>(pca_result);
    RealMatrix U = boost::get<1>(pca_result);

    if (length_normalize)
      for (uint j=0; j<s.rows(); ++j)
        s[j] /= blocksize;


    powers.push_back(sa_sum / sum(s));
    coverlaps.push_back(covarianceOverlap(sa, Ua, s, U));
  }

  TimeSeries<double> covt(coverlaps);
  TimeSeries<double> powt(powers);
  return( Datum(covt.average(), covt.variance(), powt.average()) );
}




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  int k=1;

  if (argc != 6) {
    cerr << "Usage- bcom model traj sel [1=local avg|0=global avg] blocks\n";
    exit(0);
  }

  AtomicGroup model = createSystem(argv[k++]);
  pTraj traj = createTrajectory(argv[k++], model);

  AtomicGroup subset = selectAtoms(model, argv[k++]);
  int local_flag = atoi(argv[k++]);
  vector<uint> blocksizes = parseRangeList<uint>(argv[k++]);

  vector<AtomicGroup> ensemble;
  readTrajectory(ensemble, subset, traj);
 
  // First, get the complete PCA result...
  boost::tuple<std::vector<XForm>, greal, int> ares = iterativeAlignment(ensemble);
  AtomicGroup avg = averageStructure(ensemble);
  NoAlignPolicy policy(avg, local_flag);
  boost::tuple<RealMatrix, RealMatrix> res = pca(ensemble, policy);

  RealMatrix Us = boost::get<0>(res);
  RealMatrix UA = boost::get<1>(res);

  if (length_normalize)
    for (uint i=0; i<Us.rows(); ++i)
      Us[i] /= traj->nframes();

  cout << "# " << hdr << endl;
  cout << "# Config flags: length_normalize=" << length_normalize << endl;
  cout << "# Alignment converged to " << boost::get<1>(ares) << " in " << boost::get<2>(ares) << " iterations\n";
  cout << "# n\tCoverlap\tVariance\tAvg Pow Ratio\n";
  // Now iterate over all requested block sizes...
  
  PercentProgress watcher;
  ProgressCounter<PercentTrigger, EstimatingCounter> slayer(PercentTrigger(0.1), EstimatingCounter(blocksizes.size()));
  slayer.attach(&watcher);
  slayer.start();

  for (vector<uint>::iterator i = blocksizes.begin(); i != blocksizes.end(); ++i) {
    Datum result = blocker(UA, Us, ensemble, *i, policy);
    cout << *i << "\t" << result.avg_coverlap << "\t" << result.var_coverlap << "\t" << result.avg_power << endl;
    slayer.update();
  }

  slayer.finish();

}
