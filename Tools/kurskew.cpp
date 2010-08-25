/*
  kurskew

  Compute skew and kurtosis for each column in a matrix
*/


/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2010 Tod D. Romo
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester

  This package (LOOS) is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation under version 3 of the License.

  This package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/



#include <loos.hpp>

using namespace std;
using namespace loos;


double mean(const DoubleMatrix& M, const uint i) {

  double m = 0.0;
  for (uint j=0; j<M.rows(); ++j)
    m += M(j, i);
  return(m / M.rows());
}


double stddev(const DoubleMatrix& M, const uint i, const double m) {

  double s = 0.0;
  for (uint j=0;j<M.rows(); ++j) {
    double d = M(j, i) - m;
    s += d*d;
  }

  s /= (M.rows() - 1);
  return(sqrt(s));
}


double moment(const DoubleMatrix& M, const uint i, const double m, const double s, const double p) {

  double x = 0.0;
  for (uint j=0; j<M.rows(); ++j) {
    double d = (M(j, i) - m) / s;
    x += pow(d, p);
  }
  return(x / M.rows());
}


int main(int argc, char *argv[]) {
  
  if (argc == 1) {
    cout << "Usage- kurskew matrix >output\n";
    exit(0);
  }

  string hdr = invocationHeader(argc, argv);

  DoubleMatrix M;
  readAsciiMatrix(argv[1], M);

  DoubleMatrix K(M.cols(), 2);

  for (uint i=0; i<M.cols(); ++i) {
    double avg = mean(M, i);
    double dev = stddev(M, i, avg);
    
    K(i, 0) = moment(M, i, avg, dev, 3);
    K(i, 1) = moment(M, i, avg, dev, 4) - 3.0;
  }

  writeAsciiMatrix(cout, K, hdr);
}
