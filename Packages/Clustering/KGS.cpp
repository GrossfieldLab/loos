#include "Clustering.hpp"
#include <iostream>

using std::vector;
using std::cout;
using std::endl;

using namespace Eigen;

namespace Clustering
{
  // call this to search for a cutoff stage in clustering.
  uint KGS::cutoff()
  {
    double min = avgSpread.minCoeff();
    double max = avgSpread.maxCoeff();
    double norm = (eltCount - 2) / (max - min);
#ifdef DEBUG
    cout << "penalties:" << endl
        << penalties << endl;
    cout << "avgSpreads:  " << endl
        << avgSpread << endl;
#endif
    VectorXd normAvSp = (norm * (avgSpread.array() - min) + 1).matrix();
#ifdef DEBUG
    cout << "normalized avgSpreads:" << endl
        << normAvSp << endl;
#endif
    penalties += normAvSp;
#ifdef DEBUG
    cout << "penalties after adding normAvSpreads:" << endl
        << penalties << endl;
#endif
    uint minIndex;
    penalties.minCoeff(&minIndex);
    // need to increment minindex to correspond to stage,
    // since avgSpread (and therefore penalty) undefined at stage 0.
    // This has caused us to use a penalties vector that is eltCount-1 long.
    return minIndex + 1;
  }

  void KGS::penalty()
  {
    // look up merged clustersize so we can assess change in spread.
    uint sizeA = (clusterTraj[stage - 1][minRow]).size();
    uint sizeB = (clusterTraj[stage - 1][minCol]).size();
    uint sizeAB = sizeA + sizeB;
#ifdef DEBUG
    cout << "sizeA:  " << sizeA << endl;
    cout << "sizeB:  " << sizeB << endl;
#endif
    double normSpA{0};
    double normSpB{0};
    double sumCrossDists = sizeA * sizeB * distOfMerge(stage);
    // spread of A will be sum of distances of elts in a divided by (N*(N-1)/2)
    // This is for both A and B, hence two goes to the numerator of their sum.
    // nClusters goes up to record addition of one merged (nontrivial cluster)
    if (merged)
    {
      // determine if the merge created a nontrivial cluster
      if (sizeA == 1)
        currentClusterCount++;
      else
        normSpA = 0.5 * (sizeA * (sizeA - 1)) * spreads(minRow);
      // account for the case where the merged cluster was nontrivial
      if (sizeB > 1)
      {
        currentClusterCount--;
        normSpB = 0.5 * (sizeB * (sizeB - 1)) * spreads(minCol);
      }
      // spreads(minRow) = 2*(2*(normSpA + normSpB) + sumCrossDists)/(sizeAB*(sizeAB-1));
      // remove spreads[minCol]
      removeRow(spreads, minCol);
      if (minCol < minRow)
        spreads(minRow - 1) = 2 * (2 * (normSpA + normSpB) + sumCrossDists) / (sizeAB * (sizeAB - 1));
      else
        spreads(minRow) = 2 * (2 * (normSpA + normSpB) + sumCrossDists) / (sizeAB * (sizeAB - 1));
    }
    else
    {
      // determine if the merge created a nontrivial cluster
      if (sizeB == 1)
        currentClusterCount++;
      else
        normSpB = 0.5 * (sizeB * (sizeB - 1)) * spreads(minCol);
      // account for the case where the merged cluster was nontrivial
      if (sizeA > 1)
      {
        currentClusterCount--;
        normSpA = 0.5 * (sizeA * (sizeA - 1)) * spreads(minCol);
      }
      // spreads(minCol) = 2*(2*(normSpA + normSpB) + sumCrossDists)/(sizeAB*(sizeAB-1));
      // remove spreads[minRow]
      removeRow(spreads, minRow);
      if (minRow < minCol)
        spreads(minCol - 1) = 2 * (2 * (normSpA + normSpB) + sumCrossDists) / (sizeAB * (sizeAB - 1));
      else
        spreads(minCol) = 2 * (2 * (normSpA + normSpB) + sumCrossDists) / (sizeAB * (sizeAB - 1));
    }
    // from paper, divide only by number of nontrivial clusters.
    avgSpread(stage - 1) = spreads.sum() / currentClusterCount;
    // set penalties at the number of clusters, which is the same as eltCount - stage
    penalties(stage - 1) = eltCount - stage;
  }
}