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


#include "fitter.hpp"

namespace ENM {


  double Fitter::operator()(const std::vector<double>& v) {
    enm_->setParams(v);
    if (! enm_->validParams())
      return(std::numeric_limits<double>::max());
    enm_->solve();
    
    uint n = (enm_->eigenvalues()).rows();
    uint m = (enm_->eigenvectors()).rows();

    loos::DoubleMatrix s(n-6,1);
    loos::DoubleMatrix U(m, n-6);
    
    for (uint i=0; i<n-6; ++i) {
      s[i] = 1.0 / ((enm_->eigenvalues())[n-i-1]);
      for (uint j=0; j<m; ++j)
        U(j, i) = (enm_->eigenvectors())(j, n-i-1);
    }

    if (normalize_) {
      double scale = normalizePower(s);
      for (uint j=0; j<s.rows(); ++j)
        s[j] *= scale;
    }

    double d = loos::Math::covarianceOverlap(s, U, ref_eigvals_, ref_eigvecs_);

    if (verbose_) {
      std::cout << name_ << ": ";
      std::cout << "\t(";
      for (uint i=0; i<v.size(); ++i)
        std::cout << v[i] << (i == v.size()-1 ? "" : ",");
      std::cout << ") = " << d << std::endl;
    }

    // Maximizing covariance overlap, remember?
    return(-d);
  }


  double Fitter::normalizePower(const loos::DoubleMatrix& s) {
    double sumA = 0.0;
    for (uint j=0; j<s.rows(); ++j)
      sumA += s[j];

    double sumB = 0.0;
    for (uint j=0; j<ref_eigvals_.rows(); ++j)
      sumB += ref_eigvals_[j];

    return(sumB/sumA);
  }


};
