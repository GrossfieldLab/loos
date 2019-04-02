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
MatrixXd readMatrixFromStream(istream& input){
  vector<vector<double>> matbuff;
  string line;
  double elt;
  while(getline(input, line)){
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
  TriangularBase<MatrixXd> result;
  for (uint i = 0; i < matbuff.size(); i++)
    for (uint j = i; j < matbuff[0].size(); j++)
      result(i,j) = matbuff[i][j];
  
  return result;
};

// takes a nxd data matrix, returns an nxn matrix containing pairwise distances
// use formula (a - b)^2 = a^2 + b^2 -2a*b.
MatrixXd pairwise_dists(const Ref<const MatrixXd> &data)
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
PermutationMatrix<Dynamic, Dynamic> sort_permutation(const Ref<const VectorXd> &v){
  // initialize original index locations
  PermutationMatrix<Dynamic, Dynamic> p(v.size());
  p.setIdentity();
  // sort indexes based on comparing values in v
  sort(p.indices().data(), p.indices().data() + p.indices().size(),
       [&v](size_t i1, size_t i2) { return v.data()[i1] < v.data()[i2]; });
  return p;
}

#endif
