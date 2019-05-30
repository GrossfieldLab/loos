#if !defined(LOOS_CLUSTER_HPP)
#define LOOS_CLUSTER_HPP
#include <eigen3/Eigen/Dense>
#include <vector>
#include <iosfwd>
#include <memory>

namespace Clustering
{
  #include "ClusteringUtils.hpp"
  #include "HAC.hpp"
  #include "AverageLinkage.hpp"
  #include "KGS.hpp"
}
#endif
