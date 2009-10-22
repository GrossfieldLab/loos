/*
  MatrixOps.hpp
*/

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo, Alan Grossfield
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

#if !defined MATRIXOPS_HPP
#define MATRIXOPS_HPP

#include <loos_defs.hpp>
#include <MatrixImpl.hpp>
#include <stdexcept>



namespace loos {


  typedef Math::Matrix<float, Math::ColMajor> RealMatrix;
  typedef Math::Matrix<double, Math::ColMajor> DoubleMatrix;

  namespace Math {


    boost::tuple<RealMatrix, RealMatrix, RealMatrix> svd(RealMatrix& M);
    boost::tuple<DoubleMatrix, DoubleMatrix, DoubleMatrix> svd(DoubleMatrix& M);
    DoubleMatrix MMMultiply(const DoubleMatrix& A, const DoubleMatrix& B, const bool transa = false, const bool transb = false);
    DoubleMatrix invert(DoubleMatrix& A, const double eps = 1e-6);

    void operator+=(DoubleMatrix& A, const DoubleMatrix& B);
    DoubleMatrix operator+(const DoubleMatrix& A, const DoubleMatrix& B);
    void operator-=(DoubleMatrix& A, const DoubleMatrix& B);
    DoubleMatrix operator-(const DoubleMatrix& A, const DoubleMatrix& B);
    void operator*=(DoubleMatrix& A, const double d);
    DoubleMatrix operator*(const DoubleMatrix& A, const double d);
    void operator*=(DoubleMatrix& A, const DoubleMatrix& B);
    DoubleMatrix operator*(const DoubleMatrix& A, const DoubleMatrix& B);
    DoubleMatrix operator-(DoubleMatrix& A);
    DoubleMatrix eye(const uint n);


  };


};


#endif
