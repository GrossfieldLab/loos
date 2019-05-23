#if !defined(LOOS_CLUSTER_HPP)
#define LOOS_CLUSTER_HPP

#include <boost/program_options.hpp>
#include <boost/format.hpp>
// #include <loos.hpp>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <fstream>
#include "kgs.hpp"

// namespace opts = loos::OptionsFramework;
// namespace po = boost::program_options;


// takes an istream containing an ascii matrix,
// returns arb. dimension matrix containing its contents
// Note: assumes matrix is triangular (since similarity scores
// for clustering must be reflexive...)
namespace clustering
{
  Eigen::MatrixXd readMatrixFromStream(std::istream &input, char commentChar = '#');

  // takes a nxd data matrix (where d is the dimensionality of the data),
  // returns an nxn matrix containing pairwise distances
  Eigen::MatrixXd pairwiseDists(const Eigen::Ref<const Eigen::MatrixXd> &data);


  // for exemplars defined as having the minimum average distance within cluster
  // Takes a vector of vectors of uints which are the cluster indexes, and a corresponding (full) distance matrix
  // Returns a vector of indexes to the minimum average distance element from each cluster. 
  std::vector<uint> getExemplars(std::vector<std::vector<uint>> &clusters, const Eigen::Ref<const Eigen::MatrixXd> &distances);

  // from <https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-     track-of-indexes>
  // provides a sort index in ASCENDING order. Apply using matrix product
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> sort_permutation(const Eigen::Ref<const Eigen::VectorXd> &v);

  // helper functions for adding and subtracting rows. Can GO AWAY with eigen3.4.
  // as of 4/2/19 that's months away, though the feature is finished and in devel.
  template <typename Derived>
  void removeRow(Eigen::PlainObjectBase<Derived> &matrix, unsigned int rowToRemove)
  {
    unsigned int numRows = matrix.rows() - 1;
    unsigned int numCols = matrix.cols();

    if (rowToRemove < numRows)
      matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.bottomRows(numRows - rowToRemove);

    matrix.conservativeResize(numRows, numCols);
  }

  template <typename Derived>
  void removeCol(Eigen::PlainObjectBase<Derived> &matrix, unsigned int colToRemove);

  // Abstract class for hierarchical agglomerative clustering.
  // Specific comparison methods inherit from here.
  class HAC
  {
    // no private data, since this class exists to provide inheritance.
  public:
    HAC(const Eigen::Ref<Eigen::MatrixXd> &e) : clusterDists(e.selfadjointView<Eigen::Upper>()),
                                  eltCount{e.cols()},
                                  distOfMerge(e.cols()) {}

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> clusterDists;
    // record a trajectory of the clustering so that you can write dendrograms or similar if desired.
    // These will all be of length matching clustering steps (Nelts-1)
    Eigen::VectorXd distOfMerge;
    // holds total number of elements to be clustered (and thus number of steps)
    uint eltCount;

    // These members change each step.
    // these will store the indexes of the coefficients sought.
    uint minRow, minCol, stage;
    // this bool stores outcome of 'merge'
    bool merged;
    // track the 'trajectory' of the clustering process
    std::vector<std::vector<std::vector<uint>>> clusterTraj;
    // the vector of pointers to each cluster at the current stage.
    // each element of cluster list will be currStg at stage == index.
    std::vector<std::unique_ptr<std::vector<uint>>> currStg;

    // need to fill this in for each type of
    virtual Eigen::RowVectorXd dist(uint A, uint B) {}
    // define a penalty function to score each level of the hierarchy.
    virtual void penalty() {}

    // Merge two clusters into whichever is larger.
    // Return true if new composite cluster is minRow, else return false
    // In the case where clusters are of equal size, merge into minRow.
    virtual bool merge();
    

    // Run through the clustering cycle, populating the 'trajectory' vectors.
    void cluster();
    // write clusters as YAML for easy later interp.
    void writeClusters(uint optStg, std::ostream &out);
  };

  // average linkage class for hierarchical clustering.
  // derive specific examples of average linkage HAC from here.
  // By definition they should all need this distance function.
  class AverageLinkage : public HAC
  {
  public:
    AverageLinkage(const Eigen::Ref<Eigen::MatrixXd> &e) : HAC(e) {}
    // this should be a terminal definition
    virtual Eigen::RowVectorXd dist(uint idxA, uint idxB);
  };
}
#endif
