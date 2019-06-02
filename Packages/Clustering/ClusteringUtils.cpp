#include "ClusteringUtils.hpp"
#include "ClusteringTypedefs.hpp"
#include <iostream>
#include <algorithm>
#include <fstream>

using std::cout;
using std::endl;
using std::istream;
using std::ostream;
using std::sort;
using std::string;
using std::stringstream;
using std::vector;

using namespace Eigen;
namespace Clustering
{
// MatrixXd
// readMatrixFromStream(istream &input, const char commentChar)
// {
//   vector<vector<double>> matbuff;
//   string line;
//   double elt;
//   while (getline(input, line))
//   {
//     // skip commets. Only permits comments at the beginning of lines.
//     if (line[0] == commentChar)
//       continue;
//     stringstream streamline(line);
//     vector<double> row;
//     // process a row here. Should work for whitespace delimited...
//     while (streamline >> elt)
//       // if a single line comment char is found, break out to line loop
//       row.push_back(elt);
//     // push the vector into the matrix buffer.
//     matbuff.push_back(row);
//   }

//   // Populate matrix with numbers.
//   // should be a better way to do this with Eigen::Map...
//   // though nb mapped eigen matricies are not the same as eigen dense mats.
//   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
//       result(matbuff[0].size(), matbuff.size());
//   for (idxT i = 0; i < matbuff.size(); i++)
//     for (idxT j = i; j < matbuff[0].size(); j++)
//       result(i, j) = matbuff[i][j];

//   return result;
// };

// use formula (a - b)^2 = a^2 + b^2 -2a*b.
// MatrixXd
// pairwiseDists(const Ref<const MatrixXd> &data)
// {
//   const VectorXd data_sq = data.rowwise().squaredNorm();
//   MatrixXd distances;
//   distances = data_sq.rowwise().replicate(data.rows()) +
//               data_sq.transpose().colwise().replicate(data.rows()) -
//               2. * data * data.transpose();
//   distances.diagonal().setZero(); // prevents nans from occurring along diag.
//   distances = distances.cwiseSqrt();
//   return distances;
// }

// possibly naive implementation, relies on keeping full similarity matrix.
// vector<idxT>
// getExemplars(vector<vector<idxT>> &clusters,
//              const Ref<const MatrixXd> &distances)
// {
//   vector<idxT> exemplars(clusters.size());
//   for (idxT cdx = 0; cdx < clusters.size(); cdx++)
//   {
//     MatrixXd clusterDists(clusters[cdx].size(), clusters[cdx].size());
//     for (idxT i = 0; i < clusters[cdx].size(); i++)
//     {
//       for (idxT j = 0; j < i; j++)
//       {
//         clusterDists(i, j) = distances(clusters[cdx][i], clusters[cdx][j]);
//       }
//     }
//     idxT centeridx;
//     clusterDists = clusterDists.selfadjointView<Upper>();
//     clusterDists.colwise().mean().minCoeff(&centeridx);
//     exemplars[cdx] = clusters[cdx][centeridx];
//   }
//   return exemplars;
// }
// from
// <https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes>
// PermutationMatrix<Dynamic, Dynamic>
// sort_permutation(const Ref<const VectorXd> &v)
// {
//   // initialize original index locations
//   PermutationMatrix<Dynamic, Dynamic> p(v.size());
//   p.setIdentity();
//   // sort indexes based on comparing values in v
//   sort(p.indices().data(),
//        p.indices().data() + p.indices().size(),
//        [&v](size_t i1, size_t i2) { return v.data()[i1] < v.data()[i2]; });
//   return p;
// }

} // namespace Clustering