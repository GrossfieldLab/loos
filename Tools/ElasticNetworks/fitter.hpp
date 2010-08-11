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

//! Class for fitting ENM spring parameters by comparing an ENM and PCA results
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


  double operator()(const std::vector<double>& v);

private:

  // Scale factor to make the power contained in eigenvalues s match
  // the reference eigenvalues.  Assumes zero eigenpairs have already
  // been trimmed off.

  double normalizePower(const loos::DoubleMatrix& s);

private:
  ElasticNetworkModel *enm_;
  loos::DoubleMatrix ref_eigvals_;
  loos::DoubleMatrix ref_eigvecs_;

  bool normalize_;
  bool verbose_;
  std::string name_;


};




//! Combines multiple ENMFitters together to return a joint overlap
class FitAggregator {
public:
  FitAggregator() : iters_(0), verbose_(true) { }

  bool verbose() const { return(verbose_); }
  void verbose(const bool b) { verbose_ = b; }

  uint iterations() const { return(iters_); }


  void push_back(ENMFitter* p) { fitters.push_back(p); }

  double operator()(const std::vector<double>& v) {
    double sum = 0.0;

    for (std::vector<ENMFitter*>::iterator i = fitters.begin(); i != fitters.end(); ++i)
      sum += (**i)(v);

    sum /= fitters.size();
    
    ++iters_;
    if (verbose_) 
      std::cout << "* (" << iters_ << ") Joint = " << -sum << "\n";

    return(sum);
  }

  void resetCount() { iters_ = 0; }


private:
  uint iters_;
  bool verbose_;
  std::vector<ENMFitter*> fitters;

};



#endif
