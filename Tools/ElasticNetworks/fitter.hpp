/*
 Fitter for ENM parameters
 (c) 2010 Tod D. Romo, Grossfield Lab, URMC


 Notes:

 It's assumed that we will always be comparing with PCA results and
 the number of skipped eigenpairs (for both sides) is 6.


*/



#if !defined(FITTER_HPP)
#define FITTER_HPP

#include <loos.hpp>
#include "enm-lib.hpp"


class ENMFitter {
public:
  ENMFitter(ElasticNetworkModel* model, const loos::DoubleMatrix& s, const loos::DoubleMatrix& U) :
    enm_(model),
    normalize_(false)
  {
    uint m = U.rows();
    uint n = s.rows();

    ref_eigvals_ = submatrix(s, Range(0, n-6), Range(0, 1));
    ref_eigvecs_ = submatrix(U, Range(0, m), Range(0, n-6));


    // These PCA eigenpairs actually come from an SVD, so must square
    // the svals to make eigenvalues out of them...

    for (uint j=0; j<n-6; ++j)
      ref_eigvals_[j] *= ref_eigvals_[j];
  }

  void normalize(const bool b) { normalize_ = b; }
  bool normalize() const { return(normalize_); }


  double operator()(const std::vector<double>& v) {
    enm_->setConstants(v);
    enm_->solve();
    
    uint n = (enm_->eigenvalues()).rows();
    uint m = (enm_->eigenvectors()).rows();

    loos::DoubleMatrix s = submatrix(enm_->eigenvalues(), Range(6, n), Range(0, 1));
    loos::DoubleMatrix U = submatrix(enm_->eigenvectors(), Range(0, m), Range(6, n));

    if (normalize_) {
      double scale = normalizePower(s);
      for (uint j=0; j<s.rows(); ++j)
        s[j] *= scale;
    }

    return(loos::Math::covarianceOverlap(s, U, ref_eigvals_, ref_eigvecs_));
  }


private:

  // Scale factor to make the power contained in eigenvalues s match
  // the reference eigenvalues.  Assumes zero eigenpairs have already
  // been trimmed off.

  double normalizePower(const loos::DoubleMatrix& s) {
    double sumA = 0.0;
    for (uint j=0; j<s.rows(); ++j)
      sumA += s[j];

    double sumB = 0.0;
    for (uint j=0; j<ref_eigvals_.rows(); ++j)
      sumB += ref_eigvals_[j];

    return(sumB/sumA);
  }

private:
  ElasticNetworkModel *enm_;
  loos::DoubleMatrix ref_eigvals_;
  loos::DoubleMatrix ref_eigvecs_;

  bool normalize_;


};



#endif
