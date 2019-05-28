#ifndef LOOS_HAC_HPP
#define LOOS_HAC_HPP
#include <eigen3/Eigen/Dense>

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
  virtual Eigen::RowVectorXd dist(uint A, uint B);
  // define a penalty function to score each level of the hierarchy.
  virtual void penalty();

  // Merge two clusters into whichever is larger.
  // Return true if new composite cluster is minRow, else return false
  // In the case where clusters are of equal size, merge into minRow.
  virtual bool merge();


  // Run through the clustering cycle, populating the 'trajectory' vectors.
  void cluster();

};
#endif
