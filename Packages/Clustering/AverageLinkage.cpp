#include "AverageLinkage.hpp"

namespace Clustering 
{
  Eigen::RowVectorXd AverageLinkage::dist(idxT idxA, idxT idxB)
  {
    idxT sizeA = currStg[idxA]->size();
    idxT sizeB = currStg[idxB]->size();
    return (sizeA * clusterDists.row(idxA) + sizeB * clusterDists.row(idxB)) / (sizeA + sizeB);
  }
}