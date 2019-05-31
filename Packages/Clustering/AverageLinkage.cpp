#include "AverageLinkage.hpp"

namespace Clustering 
{
  Eigen::RowVectorXd AverageLinkage::dist(uint idxA, uint idxB)
  {
    uint sizeA = currStg[idxA]->size();
    uint sizeB = currStg[idxB]->size();
    return (sizeA * clusterDists.row(idxA) + sizeB * clusterDists.row(idxB)) / (sizeA + sizeB);
  }
}