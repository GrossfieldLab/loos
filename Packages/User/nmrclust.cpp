// nmrclust.cpp follows Kelly, Gardner, and Sutcliffe, Prot. Eng. 9 11 1063-1065 (1996) 
// hereafter KGS
// To perform exactly the analysis specified there, one must first apply
// one of the all-to-all rmsd tolls (such as rmsds or multi-rmsds) before
// running this script. Those tools write their RMSD matricies to cout, and
// this script reads from cin, so the effect can be achieved through a pipe.
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "cluster.hpp"
using namespace Eigen;
using namespace std;

class NMRClust: public AverageLinkage {
public:
  NMRClust(const Ref<MatrixXd> &e) : AverageLinkage(e),
                                     avgSpread(e.rows()-1),
                                     spreads(e.rows()),
                                     nClusters(e.rows()-1),
                                     penalties(e.rows()-1) {}
  // need to track the average spread at each stage of the clustering.  
  VectorXd avgSpread;
  // need to track the number of nontrivial clusters at each stage
  VectorXi nClusters;
  // compute penalties for each step
  VectorXd penalties;

  // call this to search for a cutoff stage in clustering.
  uint cutoff()
  { 
    double min = avgSpread.minCoeff();
    double norm = (avgSpread.size()-1)/(avgSpread.maxCoeff()-min);
    penalties = norm*(avgSpread.array() - min) + 1;
    uint minIndex;
    penalties.minCoeff(&minIndex);
    return minIndex;
  }

private:
  // this will change per round of clustering
  VectorXd spreads; 
  void penalty()
  { 
    double currentClusterCount = nClusters[stage-1];
    uint sizeA = clusterList[minRow]->size();
    uint sizeB = clusterList[minCol]->size();
    double sumCrossDists = sizeA*sizeB*distOfMerge[stage];
    double newClusterSpread = 2*(spreads[minRow]/(sizeA*(sizeA-1)) + spreads[minCol]/(sizeB*(sizeB-1))) + sumCrossDists;
    // nClusters goes up to record addition of one merged (nontrivial cluster)
    currentClusterCount ++;
    if (merged)
    { // accout for the case where the merged cluster was also nontrivial
      if (sizeB > 1)
        currentClusterCount --;
      spreads[minRow] = newClusterSpread;
      // remove spreads[minCol]
      spreads[minCol] = 0;
    }
    else
    { // account for the case where the merged cluster was also nontrivial
      if (sizeA > 1)
        currentClusterCount --; 
      spreads[minCol] = newClusterSpread;
      // remove spreads[minRow]
      spreads[minRow] = 0;
    }
    nClusters[stage] = currentClusterCount;
    avgSpread[stage] = spreads.sum()/currentClusterCount;
  }
};

int main()
{
  MatrixXd similarity_scores;
  // read similarity scores into matrix, this is step 1 in CGS fig. 1
  
}