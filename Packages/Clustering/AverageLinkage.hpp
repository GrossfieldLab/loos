#ifndef LOOS_AVG_LINK_HPP
#define LOOS_AVG_LINK_HPP
#include "ClusteringTypedefs.hpp"
#include "HAC.hpp"

// average linkage class for hierarchical clustering.
// derive specific examples of average linkage HAC from here.
// By definition they should all need this distance function.
namespace Clustering
{
class AverageLinkage : public HAC
{
public:
  AverageLinkage(const Eigen::Ref<Eigen::Matrix<dtype, Eigen::Dynamic, Eigen::Dynamic>> &e) : HAC(e) {}
  // this should be a terminal definition
  Eigen::Matrix<dtype, 1, Eigen::Dynamic> dist(idxT idxA, idxT idxB);
};
} // namespace Clustering
#endif