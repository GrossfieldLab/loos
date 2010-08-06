//
// vsa-lib
//
// (c) 2010 Tod D. Romo, Grossfield Lab, URMC
//


#if !defined(VSA_LIB_HPP)
#define VSA_LIB_HPP

#include <loos.hpp>
#include "enm-lib.hpp"

#if defined(__linux__)
extern "C" {
  void dsygvx_(int*, char*, char*, char*, int*, double*, int*, double*, int*, double*, double*, int*, int*, double*, int*, double*, double*, int*, double*, int*, int*, int*, int*);
  void dpotrf_(char*, int*, double*, int*, int*);
  void dtrmm_(char*, char*, char*, char*, int*, int*, double*, double*, int*, double*, int*);
}
#endif




class VSA : public ElasticNetworkModel {
public:
  VSA(SuperBlock* blocker, const uint subn) : 
    ElasticNetworkModel(blocker),
    subset_size_(subn)
  {}

  VSA(SuperBlock* blocker, const uint subn, const loos::DoubleMatrix& M) :
    ElasticNetworkModel(blocker),
    subset_size_(subn),
    masses_(M)
  {}



  void solve();
  

  void setMasses(const loos::DoubleMatrix& M) {
    masses_ = M;
  };


  // Free up internal storage...
  void free() {
    masses_.reset();
    Msp_.reset();
    Hssp_.reset();
  }


private:
  boost::tuple<loos::DoubleMatrix, loos::DoubleMatrix> eigenDecomp(loos::DoubleMatrix& A, loos::DoubleMatrix& B);
  loos::DoubleMatrix massWeight(loos::DoubleMatrix& U, loos::DoubleMatrix& M);

  
private:
  uint subset_size_;
  loos::DoubleMatrix masses_;

  loos::DoubleMatrix Msp_;
  loos::DoubleMatrix Hssp_;
};


#endif
