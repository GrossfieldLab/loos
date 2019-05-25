#ifndef AVG_LINK
#define AVG_LINK

#include "Clustering.hpp"

// average linkage class for hierarchical clustering.
// derive specific examples of average linkage HAC from here.
// By definition they should all need this distance function.
class AverageLinkage : public Clustering::HAC
{
public:
  AverageLinkage(const Eigen::Ref<Eigen::MatrixXd> &e) : Clustering::HAC(e) {}
  // this should be a terminal definition
  Eigen::RowVectorXd dist(uint idxA, uint idxB);
};
#endif