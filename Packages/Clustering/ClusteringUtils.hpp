// ClusteringUtils.hpp
#ifndef LOOS_CLUSTERING_UTILS
#define LOOS_CLUSTERING_UTILS
#include <Eigen/Dense>
#include <iosfwd>

// takes an istream containing an ascii matrix,
// returns arb. dimension matrix containing its contents
// Note: assumes matrix is triangular (since similarity scores
// for clustering must be reflexive...)
Eigen::MatrixXd readMatrixFromStream(std::istream &input,
                                     char commentChar = '#');

// takes a nxd data matrix (where d is the dimensionality of the data),
// returns an nxn matrix containing pairwise distances
Eigen::MatrixXd pairwiseDists(const Eigen::Ref<const Eigen::MatrixXd> &data);

// for exemplars defined as having the minimum average distance within cluster
// Takes a vector of vectors of uints which are the cluster indexes, and a
// corresponding (full) distance matrix Returns a vector of indexes to the
// minimum average distance element from each cluster.
std::vector<uint>
getExemplars(std::vector<std::vector<uint>> &clusters,
             const Eigen::Ref<const Eigen::MatrixXd> &distances);

// write clusters as JSON for easy transport to analysis context.
void vectorVectorsAsJSONArr(std::vector<std::vector<uint>> &clusters,
                            std::ostream &out, const std::string &indent = "  ",
                            const std::string &offset = "  ");

// write single iterable (eg exemplars) as JSON for easy transport to analysis
// context.
template <typename ForLoopable>
void containerAsJSONArr(ForLoopable &container, std::ostream &out,
                        const std::string &indent = "  ",
                        const std::string &offset = "  ");

// write iterable containter to JSON arr on one line.
template <typename ForLoopable>
void containerAsOneLineJSONArr(ForLoopable &container, ostream &out);

// provides a sort index in ASCENDING order. Apply using matrix product
Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>
sort_permutation(const Eigen::Ref<const Eigen::VectorXd> &v);

// helper functions for adding and subtracting rows. Can GO AWAY with eigen3.4.
// as of 4/2/19 that's months away, though the feature is finished and in devel.
template <typename Derived>
void removeRow(Eigen::PlainObjectBase<Derived> &matrix,
               unsigned int rowToRemove);

template <typename Derived>
void removeCol(Eigen::PlainObjectBase<Derived> &matrix,
               unsigned int colToRemove);
#endif
