#ifndef LOOS_KGS_HPP
#define LOOS_KGS_HPP
#include "Clustering.hpp"

class KGS : public Clustering::AverageLinkage
{
public:
  KGS(const Eigen::Ref<MatrixXd> &e) : AverageLinkage(e),
                                       refDists(e.selfadjointView<Eigen::Upper>()),
                                       penalties(e.rows() - 1),
                                       avgSpread(e.rows() - 1),
                                       currentClusterCount{0} {}
  // need to track the average spread at each stage of the clustering.
  Eigen::VectorXd avgSpread; 
  // need to track the number of NONTRIVIAL clusters at each stage
  // => nClusters != currStg.size() except in cases where all
  // clusters are composite, which is not guaranteed until the very last stage.
  uint currentClusterCount;
  // compute penalties for each step
  Eigen::VectorXd penalties; // = VectorXd::Zero(eltCount-1);

  // Reference dists needed to back out cluster exemplars
  Eigen::MatrixXd refDists;

  // call this to search for a cutoff stage in clustering.
  uint cutoff();
  
private:
  // this will change per round of clustering
  Eigen::VectorXd spreads;
  void penalty();
};
#endif