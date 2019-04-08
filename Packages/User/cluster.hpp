#if !defined(LOOS_CLUSTER_HPP)
#define LOOS_CLUSTER_HPP

#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <loos.hpp>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <fstream>

namespace opts = loos::OptionsFramework;
namespace po = boost::program_options;

using namespace Eigen;
using namespace std;

// takes an istream containing an ascii matrix,
// returns arb. dimension matrix containing its contents
// Note: assumes matrix is triangular (since similarity scores
// for clustering must be reflexive...)
MatrixXd readMatrixFromStream(istream &input)
{
  vector<vector<double>> matbuff;
  string line;
  double elt;
  while (getline(input, line))
  {
    stringstream streamline(line);
    vector<double> row;
    // process a row here. Should work for whitespace delimited...
    while (streamline >> elt)
      row.push_back(elt);
    // push the vector into the matrix buffer.
    matbuff.push_back(row);
  }

  // Populate matrix with numbers.
  // should be a better way to do this with Eigen::Map...
  // though nb mapped eigen matricies are not the same as eigen dense mats.
  Matrix<double, Dynamic, Dynamic, RowMajor> result(matbuff[0].size(), matbuff.size());
  for (uint i = 0; i < matbuff.size(); i++)
    for (uint j = i; j < matbuff[0].size(); j++)
      result(i, j) = matbuff[i][j];

  return result;
};

// takes a nxd data matrix (where d is the dimensionality of the data), 
// returns an nxn matrix containing pairwise distances
// use formula (a - b)^2 = a^2 + b^2 -2a*b.
MatrixXd pairwiseDists(const Ref<const MatrixXd> &data)
{
  const VectorXd data_sq = data.rowwise().squaredNorm();
  MatrixXd distances;
  distances = data_sq.rowwise().replicate(data.rows()) + data_sq.transpose().colwise().replicate(data.rows()) - 2. * data * data.transpose();
  distances.diagonal().setZero(); // prevents nans from occurring along diag.
  distances = distances.cwiseSqrt();
  return distances;
}

// from <https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-     track-of-indexes>
// provides a sort index in ASCENDING order. Apply using matrix product
PermutationMatrix<Dynamic, Dynamic> sort_permutation(const Ref<const VectorXd> &v)
{
  // initialize original index locations
  PermutationMatrix<Dynamic, Dynamic> p(v.size());
  p.setIdentity();
  // sort indexes based on comparing values in v
  sort(p.indices().data(), p.indices().data() + p.indices().size(),
       [&v](size_t i1, size_t i2) { return v.data()[i1] < v.data()[i2]; });
  return p;
}

// helper functions for adding and subtracting rows. Can GO AWAY with eigen3.4.
// as of 4/2/19 that's months away, though the feature is finished and in devel.
template <typename Derived>
void removeRow(EigenBase<Derived> &matrix, unsigned int rowToRemove)
{
  unsigned int numRows = matrix.rows() - 1;
  unsigned int numCols = matrix.cols();

  if (rowToRemove < numRows)
    matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.bottomRows(numRows - rowToRemove);

  matrix.conservativeResize(numRows, numCols);
}

template <typename Derived>
void removeCol(EigenBase<Derived> &matrix, unsigned int colToRemove)
{
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols() - 1;

  if (colToRemove < numCols)
    matrix.block(0, colToRemove, numRows, numCols - colToRemove) = matrix.rightCols(numCols - colToRemove);

  matrix.conservativeResize(numRows, numCols);
}

// Abstract class for hierarchical agglomerative clustering.
// Specific comparison methods inherit from here.
class HAC {
  // no private data, since this class exists to provide inheritance.
public:
  HAC(const Ref<MatrixXd> &e) : clusterDists(e.selfadjointView<Upper>()),
                                eltCount{e.cols()},
                                distOfMerge(e.cols() - 1),
                                mergeTraj(e.cols() - 1, 2) {}


  Matrix<double, Dynamic, Dynamic, RowMajor> clusterDists;
  // record a trajectory of the clustering so that you can write dendrograms or similar if desired.
  // each merge event is a pair of cluster indices for the clusters at the stage recorded by the primary index.
  Matrix<uint, Dynamic, 2, RowMajor> mergeTraj;
  // These will all be of length matching clustering steps (Nelts-1)
  VectorXd distOfMerge;
  // holds total number of elements to be clustered (and thus number of steps)
  uint eltCount;
    
  // These members change each step.
  // these will store the indexes of the coefficients sought.
  uint minRow, minCol, stage;
  // this bool stores outcome of 'merge'
  bool merged;
  // the vector of pointers to each cluster at the current stage.
  vector<unique_ptr<vector<uint>>> clusterList;
  
  // need to fill this in for each type of 
  virtual VectorXd dist(uint A, uint B){}
  // define a penalty function to score each level of the hierarchy.
  virtual void penalty(){}
  
  // Merge two clusters, return true if merged cluster was first provided index, false otherwise.
  // In the case where clusters are of equal size, takes the first index provided.
  virtual bool merge(vector<unique_ptr<vector<uint>>> &clusterList, uint idxA, uint idxB)
  {
    bool ret;
    // hopefully these will be inlined/elided?!
    vector<uint> clusterA = *clusterList[idxA];
    vector<uint> clusterB = *clusterList[idxB];
    if (clusterA.size() < clusterB.size())
    {
      clusterB.insert(clusterB.end(), clusterA.begin(), clusterA.end());
      clusterList.erase(clusterList.begin() + idxA);
      ret = false;
    }
    else
    {
      clusterA.insert(clusterA.end(), clusterB.begin(), clusterB.end());
      clusterList.erase(clusterList.begin() + idxB);
      ret = true;
    }
    return ret;
  }

  // gets the cluster list. Will match cluster list at current step. 
  virtual vector<unique_ptr<vector<uint>>> getclusterList()
  {
    return clusterList;
  }


  // Run through the clustering cycle, populating the 'trajectory' vectors.
  virtual void cluster()
  {
    // initialize the list of cluster indices with one index per cluster
    for (uint i = 0; i < eltCount; i++)
    {
      unique_ptr<vector<uint>> cluster_ptr(new vector<uint>{i});
      clusterList.push_back(cluster_ptr);
    }
   for (stage = 0; stage < eltCount - 1; stage++)
    {
      // bind the minimum distance found for dendrogram construction
      distOfMerge[stage] = clusterDists.minCoeff(&minRow, &minCol);
      mergeTraj.row(stage) << minRow, minCol;
      // build merged row. Must happen before clusterList merge is performed.
      VectorXd mergedRow = dist(minRow, minCol);
      // merge the clusters into whichever of the two is larger. Erase the other.
      merged = merge(clusterList, minRow, minCol);
      // compute the penalty, if such is needed. Needs cluster merged into.
      penalty();
      // update the matrix of clusterDists
      if (merged)
      { // minRow was the cluster merged into
        // update clusterDists to zero out minCol column & row
        removeRow(clusterDists, minCol);
        removeCol(clusterDists, minCol);
        // recalculate minRow column and row
        if(minCol < minRow)
          minRow--;
        clusterDists.row(minRow) = mergedRow; 
        clusterDists.col(minRow) = mergedRow.transpose();
      }
      else
      { // minCol was the cluster merged into
        // update clusterDists to delete minRow column & row
        removeRow(clusterDists, minRow);
        removeCol(clusterDists, minRow);
        // recalculate minCol column and row
        if (minRow < minCol)
          minCol--;
        clusterDists.row(minCol) = mergedRow;
        clusterDists.col(minCol) = mergedRow.transpose();
      }

    }
  }
};

// average linkage class for hierarchical clustering.
// derive specific examples of average linkage HAC from here.
// By definition they should all need this distance function.
class AverageLinkage: HAC {
public:
  // this should be a terminal definition
  virtual VectorXd dist(uint idxA, uint idxB)
  { 
    uint sizeA = clusterList[idxA]->size();
    uint sizeB = clusterList[idxB]->size();
    return (sizeA*clusterDists.row(idxA) + sizeB*clusterDists.row(idxB))/(sizeA+sizeB);
  }
};

#endif
