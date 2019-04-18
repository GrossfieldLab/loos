#if !defined(LOOS_CLUSTER_HPP)
#define LOOS_CLUSTER_HPP

#include <boost/program_options.hpp>
#include <boost/format.hpp>
// #include <loos.hpp>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <fstream>

// namespace opts = loos::OptionsFramework;
// namespace po = boost::program_options;

using namespace Eigen;
using namespace std;

// takes an istream containing an ascii matrix,
// returns arb. dimension matrix containing its contents
// Note: assumes matrix is triangular (since similarity scores
// for clustering must be reflexive...)
MatrixXd readMatrixFromStream(istream &input, char commentChar = '#')
{
  vector<vector<double>> matbuff;
  string line;
  double elt;
  while (getline(input, line))
  {
    // skip commets. Only permits comments at the beginning of lines.
    if (line[0] == commentChar)
      continue;
    stringstream streamline(line);
    vector<double> row;
    // process a row here. Should work for whitespace delimited...
    while (streamline >> elt)
      // if a single line comment char is found, break out to line loop
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
void removeRow(PlainObjectBase<Derived> &matrix, unsigned int rowToRemove)
{
  unsigned int numRows = matrix.rows() - 1;
  unsigned int numCols = matrix.cols();

  if (rowToRemove < numRows)
    matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.bottomRows(numRows - rowToRemove);

  matrix.conservativeResize(numRows, numCols);
}

template <typename Derived>
void removeCol(PlainObjectBase<Derived> &matrix, unsigned int colToRemove)
{
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols() - 1;

  if (colToRemove < numCols)
    matrix.block(0, colToRemove, numRows, numCols - colToRemove) = matrix.rightCols(numCols - colToRemove);

  matrix.conservativeResize(numRows, numCols);
}

// Abstract class for hierarchical agglomerative clustering.
// Specific comparison methods inherit from here.
class HAC
{
  // no private data, since this class exists to provide inheritance.
public:
  HAC(const Ref<MatrixXd> &e) : clusterDists(e.selfadjointView<Upper>()),
                                eltCount{e.cols()},
                                distOfMerge(e.cols()) {}

  Matrix<double, Dynamic, Dynamic, RowMajor> clusterDists;
  // record a trajectory of the clustering so that you can write dendrograms or similar if desired.
  // These will all be of length matching clustering steps (Nelts-1)
  VectorXd distOfMerge;
  // holds total number of elements to be clustered (and thus number of steps)
  uint eltCount;

  // These members change each step.
  // these will store the indexes of the coefficients sought.
  uint minRow, minCol, stage;
  // this bool stores outcome of 'merge'
  bool merged;
  // track the 'trajectory' of the clustering process
  vector<vector<vector<uint>>> clusterTraj;
  // the vector of pointers to each cluster at the current stage.
  // each element of cluster list will be currStg at stage == index.
  vector<unique_ptr<vector<uint>>> currStg;

  // need to fill this in for each type of
  virtual RowVectorXd dist(uint A, uint B) {}
  // define a penalty function to score each level of the hierarchy.
  virtual void penalty() {}

  // Merge two clusters into whichever is larger.
  // Return true if new composite cluster is minRow, else return false
  // In the case where clusters are of equal size, merge into minRow.
  virtual bool merge()
  {
    bool ret;
    uint sizeA = currStg[minRow]->size();
    uint sizeB = currStg[minCol]->size();
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
    vector<vector<uint>> recordAtStg(currStg.size());
    for (uint i = 0; i < currStg.size(); i++)
    {
      recordAtStg[i] = *(currStg[i]);
    }
    clusterTraj.push_back(recordAtStg);
    return ret;
  }

  // Run through the clustering cycle, populating the 'trajectory' vectors.
  void cluster()
  {

    // initialize the list of cluster indices with one index per cluster
    vector<vector<uint>> recordCurrStg(eltCount);
    for (uint i = 0; i < eltCount; i++)
    {
      unique_ptr<vector<uint>> cluster_ptr(new vector<uint>{i});
      currStg.push_back(move(cluster_ptr));
      vector<uint> clusterRecord{i};
      recordCurrStg[i] = clusterRecord;
    }
    clusterTraj.push_back(recordCurrStg);

    // Get the max value to make the diagonal never the minCoeff (see distOfMerge[stage] below)
    double maxDist = clusterDists.maxCoeff() + 1;
    for (stage = 1; stage < eltCount; stage++)
    {
      cout << "stage:  " << stage << endl;
      // bind the minimum distance found for dendrogram construction
      distOfMerge(stage) = (clusterDists + maxDist * MatrixXd::Identity(
                              clusterDists.rows(), clusterDists.rows())
                            ).minCoeff(&minRow, &minCol);
      // build merged row. Must happen before clusterTraj merge is performed.
      VectorXd mergedRow = dist(minRow, minCol);
      // merge the clusters into whichever of the two is larger. Erase the other.
      merged = merge();
      cout << "clusters:" << endl;
      writeClusters(stage, cout);
      // compute the penalty, if such is needed. Needs cluster merged into.
      penalty();
      // update the matrix of clusterDists
      if (merged)
      { // minRow was the cluster merged into
        // update clusterDists to zero out minCol column & row
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

  void writeClusters(uint optStg, ostream &out)
  {
    cout << "# cluster_index elt_index1 elt_index2 ..."<< endl;
    for (uint i = 0; i < clusterTraj[optStg].size(); i++)
    {
      out << i << ' ';
      uint sizeOfCluster = (clusterTraj[optStg][i]).size();
      for (uint j = 0; j < sizeOfCluster; j++)
      {
        out << clusterTraj[optStg][i][j] << ' ';
      }
      out << endl;
    }
  }
};

// average linkage class for hierarchical clustering.
// derive specific examples of average linkage HAC from here.
// By definition they should all need this distance function.
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

#endif
