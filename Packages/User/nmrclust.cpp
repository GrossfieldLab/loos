// nmrclust.cpp 
// Kelly, Gardner, and Sutcliffe, Prot. Eng. 9 11 1063-1065 (1996) 
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
                                     penalties(e.rows()),
                                     currentClusterCount{0} {}
  // need to track the average spread at each stage of the clustering.  
  VectorXd avgSpread = VectorXd::Zero(eltCount);
  // need to track the number of NONTRIVIAL clusters at each stage
  // => nClusters != currStg.size() except in cases where all 
  // clusters are composite, which is not guaranteed until the very last stage.
  uint currentClusterCount;
  // compute penalties for each step
  VectorXd penalties = VectorXd::Zero(eltCount);

  // call this to search for a cutoff stage in clustering.
  uint cutoff()
  { 
    double min = avgSpread.minCoeff();
    double max = avgSpread.maxCoeff();
    double norm = (eltCount-2)/(max-min);
    cout << "avgSpreads:  "<< endl<<avgSpread << endl;
    penalties += (norm*(avgSpread.array() - min) + 1).matrix();
    cout << "penalties:"<< endl<< penalties << endl;
    uint minIndex;
    penalties.minCoeff(&minIndex);
    return minIndex;
  }

private:
  // this will change per round of clustering
  VectorXd spreads = VectorXd::Zero(eltCount); 
  void penalty()
  { 
    // look up merged clustersize so we can assess change in spread.
    uint sizeA = (clusterTraj[stage-1][minRow]).size();
    uint sizeB = (clusterTraj[stage-1][minCol]).size();
    uint sizeAB = sizeA+sizeB;
    cout << "sizeA:  " << sizeA << endl;
    double normSpA{0};
    double normSpB{0};
    double sumCrossDists = sizeA*sizeB*distOfMerge(stage);
    // spread of A will be sum of distances of elts in a divided by (N*(N-1)/2)
    // This is for both A and B, hence two goes to the numerator of their sum.
    // nClusters goes up to record addition of one merged (nontrivial cluster)
    if (merged)
    { 
      // determine if the merge created a nontrivial cluster
      if (sizeA == 1)
        currentClusterCount++;
      else
        normSpA = spreads(minRow)/(sizeA*(sizeA-1));
      // accout for the case where the merged cluster was nontrivial
      if (sizeB > 1)
      { 
        currentClusterCount--;
        normSpB = spreads(minCol)/(sizeB*(sizeB-1));
      }
      // spreads(minRow) = 2*(2*(normSpA + normSpB) + sumCrossDists)/(sizeAB*(sizeAB-1));
      // remove spreads[minCol]
      removeRow(spreads, minCol);
    }
    else
    { 
      // determine if the merge created a nontrivial cluster
      if (sizeB == 1)
        currentClusterCount++;
      else
        normSpB = spreads(minCol)/(sizeB*(sizeB-1));
      // account for the case where the merged cluster was nontrivial
      if (sizeA > 1)
      {
        currentClusterCount--; 
        normSpA = spreads(minCol)/(sizeA*(sizeA-1));
      }
      // spreads(minCol) = 2*(2*(normSpA + normSpB) + sumCrossDists)/(sizeAB*(sizeAB-1));
      // remove spreads[minRow]
      removeRow(spreads, minRow);
    }
    double ns = 2*(2*(normSpA + normSpB) + sumCrossDists)/(sizeAB*(sizeAB-1));
    spreads(minCol) = ns; 
    // from paper, divide only by number of nontrivial clusters.
    avgSpread(stage-1) = spreads.sum()/currentClusterCount;
    // set penalties at the number of clusters, which is the same as eltCount - stage 
    penalties(stage-1) = eltCount - stage;
  }
};

int main()
{
  MatrixXd similarityScores = readMatrixFromStream(cin);
  NMRClust clusterer(similarityScores);
  clusterer.cluster();
  uint optStg = clusterer.cutoff();
  clusterer.writeClusters(optStg++, cout);
}