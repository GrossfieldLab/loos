/*
  scrulib
  (c) 2010 Tod D. Romo, Grossfield Lab, URMC
*/


#if !defined(BCOMLIB_HPP)
#define BCOMLIB_HPP


#include <loos.hpp>

namespace Convergence {

  void subtractStructure(loos::RealMatrix& M, const loos::AtomicGroup& model);


  // Various policies that determine how blocks are aligned and averaged...

  struct AlignToPolicy {
    AlignToPolicy(const loos::AtomicGroup& targ) : target(targ), local_average(true) { }
    AlignToPolicy(const loos::AtomicGroup& targ, const bool flag) : target(targ), local_average(flag) { }

    loos::RealMatrix operator()(std::vector<loos::AtomicGroup>& ensemble) {
      for (std::vector<loos::AtomicGroup>::iterator i = ensemble.begin(); i != ensemble.end(); ++i)
        (*i).alignOnto(target);

      loos::RealMatrix M = loos::extractCoords(ensemble);
      if (local_average) {
        loos::AtomicGroup avg = loos::averageStructure(ensemble);
        subtractStructure(M, avg);
      } else
        subtractStructure(M, target);

      return(M);
    }

    loos::AtomicGroup target;
    bool local_average;
  };


  struct NoAlignPolicy {
    NoAlignPolicy() : local_average(true) { }
    NoAlignPolicy(const loos::AtomicGroup& avg_) : avg(avg_), local_average(false) { }
    NoAlignPolicy(const loos::AtomicGroup& avg_, const bool flag) : avg(avg_), local_average(flag) { }

    loos::RealMatrix operator()(std::vector<loos::AtomicGroup>& ensemble) {

      loos::RealMatrix M = loos::extractCoords(ensemble);
      if (local_average) {
        loos::AtomicGroup lavg = loos::averageStructure(ensemble);
        subtractStructure(M, lavg);
      } else
        subtractStructure(M, avg);
      return(M);
    }

    loos::AtomicGroup avg;
    bool local_average;
  };



  template<class ExtractPolicy>
  boost::tuple<loos::RealMatrix, loos::RealMatrix> pca(std::vector<loos::AtomicGroup>& ensemble, ExtractPolicy& extractor) {

    loos::RealMatrix M = extractor(ensemble);
    loos::RealMatrix C = loos::Math::MMMultiply(M, M, false, true);

    // Compute [U,D] = eig(C)
    char jobz = 'V';
    char uplo = 'L';
    f77int n = M.rows();
    f77int lda = n;
    float dummy;
    loos::RealMatrix W(n, 1);
    f77int lwork = -1;
    f77int info;
    ssyev_(&jobz, &uplo, &n, C.get(), &lda, W.get(), &dummy, &lwork, &info);
    if (info != 0)
      throw(loos::NumericalError("ssyev failed in loos::pca()", info));

   
    lwork = static_cast<f77int>(dummy);
    float *work = new float[lwork+1];

    ssyev_(&jobz, &uplo, &n, C.get(), &lda, W.get(), work, &lwork, &info);
    if (info != 0)
      throw(loos::NumericalError("ssyev failed in loos::pca()", info));
  
    reverseColumns(C);
    reverseRows(W);

    // Zap negative eigenvalues...
    for (uint j=0; j<W.rows(); ++j)
      if (W[j] < 0.0)
        W[j] = 0.0;

    boost::tuple<loos::RealMatrix, loos::RealMatrix> result(W, C);
    return(result);

  }

}

#endif
