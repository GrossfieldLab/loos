// ClusteringUtils.hpp
#ifndef LOOS_CLUSTERING_UTILS
#define LOOS_CLUSTERING_UTILS
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include <vector>
using std::cout;
using std::endl;
using std::istream;
using std::ostream;
using std::sort;
using std::string;
using std::stringstream;
using std::vector;

namespace Clustering
{
// takes an istream containing an ascii matrix,
// returns arb. dimension matrix containing its contents
// Note: assumes matrix is triangular (since similarity scores
// for clustering must be reflexive...)
Eigen::MatrixXd readMatrixFromStream(std::istream &input,
                                     const char commentChar = '#')
{
  using namespace Eigen;
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
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      result(matbuff[0].size(), matbuff.size());
  for (uint i = 0; i < matbuff.size(); i++)
    for (uint j = i; j < matbuff[0].size(); j++)
      result(i, j) = matbuff[i][j];

  return result;
}

// takes a nxd data matrix (where d is the dimensionality of the data),
// returns an nxn matrix containing pairwise distances
Eigen::MatrixXd pairwiseDists(const Eigen::Ref<const Eigen::MatrixXd> &data)
{
  using namespace Eigen;
  const VectorXd data_sq = data.rowwise().squaredNorm();
  MatrixXd distances;
  distances = data_sq.rowwise().replicate(data.rows()) +
              data_sq.transpose().colwise().replicate(data.rows()) -
              2. * data * data.transpose();
  distances.diagonal().setZero(); // prevents nans from occurring along diag.
  distances = distances.cwiseSqrt();
  return distances;
}

// provides a sort index in ASCENDING order. Apply using matrix product
Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>
sort_permutation(const Eigen::Ref<const Eigen::VectorXd> &v)
{
  using namespace Eigen;
  // initialize original index locations
  PermutationMatrix<Dynamic, Dynamic> p(v.size());
  p.setIdentity();
  // sort indexes based on comparing values in v
  sort(p.indices().data(),
       p.indices().data() + p.indices().size(),
       [&v](size_t i1, size_t i2) { return v.data()[i1] < v.data()[i2]; });
  return p;
}
// helper functions for adding and subtracting rows. Can GO AWAY with eigen3.4.
// as of 4/2/19 that's months away, though the feature is finished and in devel.
template <typename Derived>
void removeRow(Eigen::PlainObjectBase<Derived> &matrix,
               unsigned int rowToRemove)
{
  unsigned int numRows = matrix.rows() - 1;
  unsigned int numCols = matrix.cols();

  if (rowToRemove < numRows)
    matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) =
        matrix.bottomRows(numRows - rowToRemove);

  matrix.conservativeResize(numRows, numCols);
}

template <typename Derived>
void removeCol(Eigen::PlainObjectBase<Derived> &matrix,
               unsigned int colToRemove)
{
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols() - 1;

  if (colToRemove < numCols)
    matrix.block(0, colToRemove, numRows, numCols - colToRemove) =
        matrix.rightCols(numCols - colToRemove);

  matrix.conservativeResize(numRows, numCols);
}

// for exemplars defined as having the minimum average distance within cluster
// Takes a vector of vectors of uints which are the cluster indexes, and a
// corresponding (full) distance matrix Returns a vector of indexes to the
// minimum average distance element from each cluster.
std::vector<uint>
getExemplars(std::vector<std::vector<uint>> &clusters,
             const Eigen::Ref<const Eigen::MatrixXd> &distances)
{
  using namespace Eigen;
  vector<uint> exemplars(clusters.size());
  for (uint cdx = 0; cdx < clusters.size(); cdx++)
  {
    MatrixXd clusterDists(clusters[cdx].size(), clusters[cdx].size());
    for (uint i = 0; i < clusters[cdx].size(); i++)
    {
      for (uint j = 0; j < i; j++)
      {
        clusterDists(i, j) = distances(clusters[cdx][i], clusters[cdx][j]);
      }
    }
    uint centeridx;
    clusterDists = clusterDists.selfadjointView<Upper>();
    clusterDists.colwise().mean().minCoeff(&centeridx);
    exemplars[cdx] = clusters[cdx][centeridx];
  }
  return exemplars;
}
// write clusters as JSON for easy transport to analysis context.
template <typename Numeric>
void vectorVectorsAsJSONArr(std::vector<std::vector<Numeric>> &clusters,
                            std::ostream &out, const std::string &indent = "  ",
                            const std::string &offset = "  ")
{
  out << offset + "[" << endl;
  for (uint i = 0; i < clusters.size(); i++)
  {
    out << offset + indent + "[";
    uint lastIndex = (clusters[i]).size() - 1;
    for (uint j = 0; j < lastIndex; j++)
    {
      out << clusters[i][j] << ",";
    }
    out << clusters[i][lastIndex] << "]," << endl;
  }
  out << offset + "]";
}

// write single iterable (eg exemplars) as JSON for easy transport to analysis
// context.
template <typename ForLoopable>
void containerAsJSONArr(ForLoopable &container, std::ostream &out,
                        const std::string &indent = "  ",
                        const std::string &offset = "  ")
{
  out << "[" << endl;
  uint lastIndex = container.size() - 1;
  for (uint i = 0; i < lastIndex; i++)
  {
    out << offset + indent << container[i] << "," << endl;
  }
  out << offset << container[lastIndex] << "]";
}

// write iterable containter to JSON arr on one line.
template <typename ForLoopable>
void containerAsOneLineJSONArr(ForLoopable &container, std::ostream &out)
{
  out << "[";
  uint lastIndex = container.size() - 1;
  for (uint i = 0; i < lastIndex; i++)
  {
    out << container[i] << ",";
  }
  out << container[lastIndex] << "]";
}

} // namespace Clustering
#endif
