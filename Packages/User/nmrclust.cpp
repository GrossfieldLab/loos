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
                                     nClusters(e.rows()-1),
                                     penalties(e.rows()-1),
                                     currentClusterCount{0} {}
  // need to track the average spread at each stage of the clustering.  
  VectorXd avgSpread = VectorXd::Zero(clusterDists.rows()-1);
  // need to track the number of NONTRIVIAL clusters at each stage
  // => nClusters != currStg.size() except in cases where all 
  // clusters are composite, which is not guaranteed until the very last stage.
  VectorXi nClusters;
  uint currentClusterCount;
  // compute penalties for each step
  VectorXd penalties;

  // call this to search for a cutoff stage in clustering.
  uint cutoff()
  { 
    double min = avgSpread.minCoeff();
    double norm = (avgSpread.size()-1)/(avgSpread.maxCoeff()-min);
    penalties = norm*(avgSpread.array() - min) + 1;
    penalties += nClusters.cast<double>();
    uint minIndex;
    penalties.minCoeff(&minIndex);
    return minIndex;
  }

private:
  // this will change per round of clustering
  VectorXd spreads = VectorXd::Zero(clusterDists.rows()); 
  void penalty()
  { 
    // look up merged clustersize so we can assess change in spread.
    uint sizeA = (clusterTraj[stage-1][minRow]).size();
    uint sizeB = (clusterTraj[stage-1][minCol]).size();
    uint sizeAB = sizeA+sizeB;
    double normSpA;
    double normSpB;
    double sumCrossDists = sizeA*sizeB*distOfMerge[stage];
    // spread of A will be sum of distances of elts in a divided by (N*(N-1)/2)
    // This is for both A and B, hence two goes to the numerator of their sum.
    // nClusters goes up to record addition of one merged (nontrivial cluster)
    currentClusterCount ++;
    if (merged)
    { // accout for the case where the merged cluster was also nontrivial
      if (sizeB > 1)
      { 
        currentClusterCount--;
        normSpB = spreads[minCol]/(sizeB*(sizeB-1));
      }
      else
      {
        normSpB = 0;
      }
      normSpA = spreads[minRow]/(sizeA*(sizeA-1));
      spreads[minRow] = 2*(2*(normSpA + normSpB) + sumCrossDists)/(sizeAB*(sizeAB-1));
      // remove spreads[minCol]
      spreads[minCol] = 0;
    }
    else
    { // account for the case where the merged cluster was also nontrivial
      if (sizeA > 1)
      {
        currentClusterCount --; 
        normSpA = spreads[minCol]/(sizeA*(sizeA-1));
      }
      else
      {
        normSpA = 0;
      }
      spreads[minCol] = 2*(2*(normSpA + normSpB) + sumCrossDists)/(sizeAB*(sizeAB-1));
      // remove spreads[minRow]
      spreads[minRow] = 0;
    }
    nClusters[stage-1] = currentClusterCount;
    avgSpread[stage-1] = spreads.sum()/currentClusterCount;
  }
};

int main()
{
  MatrixXd similarityScores = readMatrixFromStream(cin);
  NMRClust clusterer(similarityScores);
  clusterer.cluster();
  uint optStg;
  clusterer.penalties.minCoeff(&optStg);
  cout << "penalties:" << endl<< clusterer.penalties << endl;
  cout << "avgSpread:" << endl << clusterer.avgSpread << endl;
  clusterer.writeClusters(optStg, cout);
}