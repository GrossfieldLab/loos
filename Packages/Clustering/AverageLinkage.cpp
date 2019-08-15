#include "AverageLinkage.hpp"

namespace Clustering 
{
  Eigen::Matrix<dtype, 1, Eigen::Dynamic> AverageLinkage::dist(idxT idxA, idxT idxB)
  {
    idxT sizeA = currStg[idxA]->size();
    idxT sizeB = currStg[idxB]->size();
    return (sizeA * clusterDists.row(idxA) + sizeB * clusterDists.row(idxB)) / (sizeA + sizeB);
  }
}