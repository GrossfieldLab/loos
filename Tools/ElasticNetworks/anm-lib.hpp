//
// (c) 2010 Tod D. Romo, Grossfield Lab, URMC
//


#if !defined(ANM_LIB_HPP)
#define ANM_LIB_HPP


#include "enm-lib.hpp"




class ANM : public ElasticNetworkModel {
public:

  void solve() {
    buildHessian();

    boost::tuple<loos::DoubleMatrix, loos::DoubleMatrix, loos::DoubleMatrix> result = svd(hessian_);
    eigenvecs_ = boost::get<0>(result);
    eigenvals_ = boost::get<1>(result);
    rsv_ = boost::get<2>(result);

    uint n = eigenvals_.rows();
    loos::reverseRows(eigenvals_);
    loos::reverseColumns(eigenvecs_);
    loos::reverseRows(rsv_);
  }


  loos::DoubleMatrix inverseHessian() {

    uint n = eigenvals_.rows();
    for (uint i=6; i<n; ++i) {
      double s = 1.0 / eigenvals_[i];
      for (uint j=0; j<n; ++j)
        rsv_(i, j) *= s;
    }

    loos::DoubleMatrix Hi = MMMultiply(rsv_, eigenvecs_, true, true);
    return(Hi);
  }


private:
  loos::DoubleMatrix rsv_;

};



#endif
