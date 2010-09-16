/*
 Fitter for ENM parameters
 (c) 2010 Tod D. Romo, Grossfield Lab, URMC


 Notes:

 It's assumed that we will always be comparing with PCA results and
 the number of skipped eigenpairs (for both sides) is 6.


*/


/** \addtogroup ENM
 *@{
 */


#if !defined(FITTER_HPP)
#define FITTER_HPP

#include <limits>
#include <loos.hpp>
#include "enm-lib.hpp"


namespace ENM {

  //! Class for fitting ENM spring parameters by comparing an ENM and PCA results
  /**
   * This class assumes that what you will be fitting are ENM results
   * against PCA results (obtained via SVD).  This means that the PCA
   * eigenvalues are expected to be singular values and must first be
   * squared.  It is also assumed that the last 6-terms are all zeros
   * (representing system rotation and translation) and will be skipped.
   *
   * Similarly, the ENM results assume that the first 6 terms will be
   * zero.
   *
   * The default behavior is to not scale the total power in the ENM to
   * match that of the PCA.
   */
  class Fitter {
  public:
    //! Associates an elastic network model with a PCA result
    /**
     * \arg \c model Pointer to an ElasticNetwork Model to fit
     * \arg \c s Single column matrix of singular values
     * \arg \c U Column-vector matrix of left singular vectors
     */
    Fitter(ElasticNetworkModel* model, const loos::DoubleMatrix& s, const loos::DoubleMatrix& U) :
      enm_(model),
      normalize_(false),
      verbose_(false)
    {
      uint m = U.rows();
      uint n = s.rows();

      ref_eigvals_ = submatrix(s, loos::Math::Range(0, n-6), loos::Math::Range(0, 1));
      ref_eigvecs_ = submatrix(U, loos::Math::Range(0, m), loos::Math::Range(0, n-6));


      // These PCA eigenpairs actually come from an SVD, so must square
      // the svals to make eigenvalues out of them...

      for (uint j=0; j<n-6; ++j)
        ref_eigvals_[j] *= ref_eigvals_[j];
    }

    //! Controls whether total power in ENM is scaled to match the PCA
    void normalize(const bool b) { normalize_ = b; }
    bool normalize() const { return(normalize_); }

    //! Name tag associated with this fit (for logging)
    void name(const std::string& s) { name_ = s; }
    std::string name() const { return(name_); }

    //! How wordy our output is
    void verbose(const bool b) { verbose_ = b; }
    bool verbose() const { return(verbose_); }


    //! Computes the covariance overlap between the ENM and the PCA
    /**
     * Takes a vector \a v of parameters to pass along to the contained
     * spring constants, then computes the ENM.  If normalization is
     * turned on, then the ENM eigenvalues are scaled so that the total
     * power is the same as the PCA.  The covariance overlap is then
     * computed and returned
     */
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




  //! Combines multiple Fitters together to return a joint (average) overlap
  class FitAggregator {
  public:
    FitAggregator() : iters_(0), verbose_(true) { }

    bool verbose() const { return(verbose_); }
    //! Determines whether or not the joint overlap is logged
    void verbose(const bool b) { verbose_ = b; }

    //! Number of total times this object has been called
    uint iterations() const { return(iters_); }

    //! Adds another system/model to fit
    void push_back(Fitter* p) { fitters.push_back(p); }

    //! Computes the joint overlap (see Fitter::operator()(const std::vector<double>& v))
    double operator()(const std::vector<double>& v) {
      double sum = 0.0;

      for (std::vector<Fitter*>::iterator i = fitters.begin(); i != fitters.end(); ++i)
        sum += (**i)(v);

      sum /= fitters.size();
    
      ++iters_;
      if (verbose_) 
        std::cout << "* (" << iters_ << ") Joint = " << -sum << "\n";

      return(sum);
    }

    //! Reset the internal call-count
    void resetCount() { iters_ = 0; }


  private:
    uint iters_;
    bool verbose_;
    std::vector<Fitter*> fitters;

  };


};


#endif



/** @} */
