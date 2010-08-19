/*
  fitter
  (c) 2010 Tod D. Romo, Grossfield Lab, URMC
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
