#include "Clustering.hpp"

using namespace Clustering;

virtual RowVectorXd dist(uint idxA, uint idxB)
{
  uint sizeA = currStg[idxA]->size();
  uint sizeB = currStg[idxB]->size();
  return (sizeA * clusterDists.row(idxA) + sizeB * clusterDists.row(idxB)) / (sizeA + sizeB);
}