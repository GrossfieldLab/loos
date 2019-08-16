#ifndef LOOS_KGS_HPP
#define LOOS_KGS_HPP
#include "ClusteringTypedefs.hpp"
#include "AverageLinkage.hpp"

namespace Clustering
{
class KGS : public AverageLinkage
{
public:
  KGS(const Eigen::Ref<Eigen::Matrix<dtype, Eigen::Dynamic, Eigen::Dynamic>> &e) : AverageLinkage(e),
                                                                                   penalties(e.rows() - 1),
                                                                                   avgSpread(e.rows() - 1),
                                                                                   currentClusterCount{0} {}

  // compute penalties for each step
  Eigen::Matrix<dtype, Eigen::Dynamic, 1> penalties; // = VectorXd::Zero(eltCount-1);

  // need to track the average spread at each stage of the clustering.
  Eigen::Matrix<dtype, Eigen::Dynamic, 1> avgSpread;

  // need to track the number of NONTRIVIAL clusters at each stage
  // => nClusters != currStg.size() except in cases where all
  // clusters are composite, which is not guaranteed until the very last stage.
  idxT currentClusterCount;

  // call this to search for a cutoff stage in clustering.
  idxT cutoff();

private:
  // this will change per round of clustering
  Eigen::Matrix<dtype, Eigen::Dynamic, 1> spreads = Eigen::Matrix<dtype, Eigen::Dynamic, 1>::Zero(eltCount);
  void penalty();
};
} // namespace Clustering
#endif