/*
 Fitter for ENM parameters
 (c) 2010 Tod D. Romo, Grossfield Lab, URMC


 Notes:

 It's assumed that we will always be comparing with PCA results and
 the number of skipped eigenpairs (for both sides) is 6.


*/



#if !defined(FITTER_HPP)
#define FITTER_HPP

#include <limits>
#include <loos.hpp>
#include "enm-lib.hpp"


class ENMFitter {
public:
  ENMFitter(ElasticNetworkModel* model, const loos::DoubleMatrix& s, const loos::DoubleMatrix& U) :
    enm_(model),
    normalize_(false),
    verbose_(false)
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

  void name(const std::string& s) { name_ = s; }
  std::string name() const { return(name_); }

  void verbose(const bool b) { verbose_ = b; }
  bool verbose() const { return(verbose_); }


  double operator()(const std::vector<double>& v) {
    if (! enm_->setConstants(v))
      return(std::numeric_limits<double>::max());
    enm_->solve();
    
    uint n = (enm_->eigenvalues()).rows();
    uint m = (enm_->eigenvectors()).rows();

    loos::DoubleMatrix s(n-6,1);
    loos::DoubleMatrix U(m, n-6);
    
    for (uint i=0; i<n-6; ++i) {
      try {
        s[i] = 1.0 / ((enm_->eigenvalues())[n-i-1]);
      }
      catch (...) {
        std::cerr << "Failure in eigenvalues at i = " << i << "\n";
        exit(-1);
      }

      for (uint j=0; j<m; ++j) {
        try {
          U(j, i) = (enm_->eigenvectors())(j, n-i-1);
        }
        catch (...) {
          std::cerr << "Failure in eigenvectors at i = " << i << ", j = " << j << "\n";
          exit(-2);
        }
      }
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
  bool verbose_;
  std::string name_;


};



#endif
