//
// Fitter for ENM parameters
// (c) 2010 Tod D. Romo, Grossfield Lab, URMC
//


#if !defined(FITTER_HPP)
#define FITTER_HPP

#include <loos.hpp>
#include "enm-lib.hpp"


class ENMFitter {
public:
  ENMFitter(const ElasticNetworkModel* model, const loos::DoubleMatrix& s, const loos::DoubleMatrix& U) :
    enm_(model),
    eigvals_(s),
    eigvecs_(U),
    normalize_(false)
  { }

  void normalize(const bool b) { normalize_ = b; }
  bool normalize() const { return(normalize_); }

  // Set the constants required for the ENM
  virtual void setConstants(const std::vector<double>& konst) =0;


private:
  ElasticNetworkModel *enm_;
  loos::DoubleMatrix eigvals_;
  loos::DoubleMatrix eigvecs_;

  bool normalize_;


};



#endif
