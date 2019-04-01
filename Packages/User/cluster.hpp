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
  MatrixXd result(matbuff.size(), matbuff[0].size());
  for (uint i = 0; i < matbuff.size(); i++)
    for (uint j = 0; j < matbuff[0].size(); j++)
      result(i,j) = matbuff[i][j];

  return result;
};

#endif
