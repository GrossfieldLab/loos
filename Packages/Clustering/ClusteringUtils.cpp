#include "Clustering.hpp"
#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>

using std::cout;
using std::ostream;
using std::istream;
using std::string;
using std::vector;
using std::endl;
using std::sort;
using std::stringstream;

using namespace Clustering;
using namespace Eigen;

MatrixXd
readMatrixFromStream(istream& input, char commentChar = '#')
{
  vector<vector<double>> matbuff;
  string line;
  double elt;
  while (getline(input, line)) {
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
};

// use formula (a - b)^2 = a^2 + b^2 -2a*b.
MatrixXd
pairwiseDists(const Ref<const MatrixXd>& data)
{
  const VectorXd data_sq = data.rowwise().squaredNorm();
  MatrixXd distances;
  distances = data_sq.rowwise().replicate(data.rows()) +
              data_sq.transpose().colwise().replicate(data.rows()) -
              2. * data * data.transpose();
  distances.diagonal().setZero(); // prevents nans from occurring along diag.
  distances = distances.cwiseSqrt();
  return distances;
}

// possibly naive implementation, relies on keeping full similarity matrix.
vector<uint>
getExemplars(vector<vector<uint>>& clusters,
             const Ref<const MatrixXd>& distances)
{
  vector<uint> exemplars(clusters.size());
  for (uint cdx = 0; cdx < clusters.size(); cdx++) {
    MatrixXd clusterDists(clusters[cdx].size(), clusters[cdx].size());
    for (uint i = 0; i < clusters[cdx].size(); i++) {
      for (uint j = 0; j < i; j++) {
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
// write clusters as a JSON array of arrays.
void
vectorVectorsAsJSONArr(vector<vector<uint>>& cluster,
                      ostream& out,
                      const string indent = "  ",
                      const string offset = "  ")
{
  out << offset + "[" << endl;
  for (uint i = 0; i < cluster.size(); i++) {
    out << offset + indent + "[";
    uint lastIndex = (cluster[i]).size() - 1;
    for (uint j = 0; j < lastIndex; j++) {
      out << cluster[i][j] << ",";
    }
    out << cluster[i][lastIndex] << "]," << endl;
  }
  out << offset + "]";
}

// write iterable container to JSON. Container must have a size() method and a
// [] operator.
template<typename ForLoopable>
void
containerAsJSONArr(ForLoopable& container,
                   ostream& out,
                   const string offset = "  ",
                   const string indent = "  ")
{
  out << "[" << endl;
  uint lastIndex = container.size() - 1;
  for (uint i = 0; i < lastIndex; i++) {
    out << offset + indent << container[i] << "," << endl;
  }
  out << offset << container[lastIndex] << "]";
}

// write iterable container to JSON arr on one line.
template<typename ForLoopable>
void
containerAsOneLineJSONArr(ForLoopable& container, ostream& out)
{
  out << "[";
  uint lastIndex = container.size() - 1;
  for (uint i = 0; i < lastIndex; i++) {
    out << container[i] << ",";
  }
  out << container[lastIndex] << "]";
}

// from
// <https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes>
PermutationMatrix<Dynamic, Dynamic>
sort_permutation(const Ref<const VectorXd>& v)
{
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
template<typename Derived>
void
removeRow(PlainObjectBase<Derived>& matrix, unsigned int rowToRemove)
{
  unsigned int numRows = matrix.rows() - 1;
  unsigned int numCols = matrix.cols();

  if (rowToRemove < numRows)
    matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) =
      matrix.bottomRows(numRows - rowToRemove);

  matrix.conservativeResize(numRows, numCols);
}

// helper functions for adding and subtracting rows. Can GO AWAY with eigen3.4.
// as of 4/2/19 that's months away, though the feature is finished and in devel.
template<typename Derived>
void
removeCol(PlainObjectBase<Derived>& matrix, unsigned int colToRemove)
{
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols() - 1;

  if (colToRemove < numCols)
    matrix.block(0, colToRemove, numRows, numCols - colToRemove) =
      matrix.rightCols(numCols - colToRemove);

  matrix.conservativeResize(numRows, numCols);
}
