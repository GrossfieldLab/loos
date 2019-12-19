#include "HAC.hpp"
#include "ClusteringUtils.hpp"
#include <iostream>

using std::cout;
using std::endl;
using std::unique_ptr;
using std::vector;
using namespace Eigen;

// Merge two clusters into whichever is larger.
// Return true if new composite cluster is minRow, else return false
// In the case where clusters are of equal size, merge into minRow.
namespace Clustering
{
bool HAC::merge()
{
  bool ret;
  idxT sizeA = currStg[minRow]->size();
  idxT sizeB = currStg[minCol]->size();
  if (sizeA < sizeB)
  {
    currStg[minCol]->insert(currStg[minCol]->end(),
                            currStg[minRow]->begin(),
                            currStg[minRow]->end());
    currStg.erase(currStg.begin() + minRow);

    ret = false;
  }
  else
  {
    currStg[minRow]->insert(currStg[minRow]->end(),
                            currStg[minCol]->begin(),
                            currStg[minCol]->end());
    currStg.erase(currStg.begin() + minCol);

    ret = true;
  }

  // append new assortment of clusters to Cluster Trajectory
  vector<vector<idxT>> recordAtStg(currStg.size());
  for (idxT i = 0; i < currStg.size(); i++)
  {
    recordAtStg[i] = *(currStg[i]);
  }
  clusterTraj.push_back(recordAtStg);
  return ret;
}

// Run through the clustering cycle, populating the 'trajectory' vectors.
void HAC::cluster()
{
  // initialize the list of cluster indices with one index per cluster
  vector<vector<idxT>> recordCurrStg(eltCount);
  for (idxT i = 0; i < eltCount; i++)
  {
    unique_ptr<vector<idxT>> cluster_ptr(new vector<idxT>{i});
    currStg.push_back(move(cluster_ptr));
    vector<idxT> clusterRecord{i};
    recordCurrStg[i] = clusterRecord;
  }
  clusterTraj.push_back(recordCurrStg);

  // Get the max value to make the diagonal never the minCoeff (see distOfMerge[stage] below)
  dtype maxDist = clusterDists.maxCoeff() + 1;
  for (stage = 1; stage < eltCount; stage++)
  {
#ifdef DEBUG
    cout << "stage:  " << stage << endl;
#endif
    // bind the minimum distance found for dendrogram construction
    distOfMerge(stage) = (clusterDists + maxDist * Matrix<dtype, Dynamic, Dynamic>::Identity(
                                                       clusterDists.rows(), clusterDists.rows()))
                             .minCoeff(&minRow, &minCol);
    // build merged row. Must happen before clusterTraj merge is performed.
    Matrix<dtype, Dynamic, 1> mergedRow = dist(minRow, minCol);
    // merge the clusters into whichever of the two is larger. Erase the other.
    merged = merge();
#ifdef DEBUG
    cout << "clusters:" << endl;
    vectorVectorsAsJSONArr(clusterTraj[stage], cout);
#endif
    // compute the penalty, if such is needed. Needs cluster merged into.
    penalty();
    // update the matrix of clusterDists
    if (merged)
    { // minRow was the cluster merged into
      // update clusterDists to zero out minCol column & row
      // removeRow<Matrix<double, -1, -1, 1>>(clusterDists, minCol);
      // removeCol<Matrix<double, -1, -1, 1>>(clusterDists, minCol);
      removeRow(clusterDists, minCol);
      removeCol(clusterDists, minCol);
      // remove the column we eliminated from our merged row of distances.
      removeRow(mergedRow, minCol);
      // recalculate minRow column and row
      if (minCol < minRow)
        minRow--;
      // Note that the dist matrix will not have a zero at this row/col after doing this.
      // because of how we have increased the values of the diagonal anyway for the mincoeff
      // this should not interfere with anything. But if you're relying on the diagonal to be zero...
      clusterDists.row(minRow) = mergedRow;
      clusterDists.col(minRow) = mergedRow.transpose();
    }
    else
    { // minCol was the cluster merged into
      // update clusterDists to delete minRow column & row
      removeRow(clusterDists, minRow);
      removeCol(clusterDists, minRow);
      // remove the column we eliminated from our merged row of distances.
      removeRow(mergedRow, minRow);
      // recalculate minCol column and row
      if (minRow < minCol)
        minCol--;

      clusterDists.row(minCol) = mergedRow;
      clusterDists.col(minCol) = mergedRow.transpose();
    }
  }
  stage--;
}
} // namespace Clustering