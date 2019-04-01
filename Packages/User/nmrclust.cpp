#include <iostream>
#include <eigen3/Eigen/Dense>
#include "cluster.hpp"
using Eigen::MatrixXd;

int main()
{
  MatrixXd m;
  m = readMatrixFromStream(std::cin);
  std::cout << m << std::endl;
}