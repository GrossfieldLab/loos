#include "Clustering.hpp"

using namespace Clustering

class AverageLinkage : public HAC
{
public:
  AverageLinkage(const Ref<MatrixXd> &e) : HAC(e) {}
  // this should be a terminal definition
  virtual RowVectorXd dist(uint idxA, uint idxB)
  {
    uint sizeA = currStg[idxA]->size();
    uint sizeB = currStg[idxB]->size();
    return (sizeA * clusterDists.row(idxA) + sizeB * clusterDists.row(idxB)) / (sizeA + sizeB);
  }
};
