/*
  bcomlib
*/




/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2010, Tod D. Romo
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

// @cond PACKAGES_INTERNAL


#if !defined(LOOS_BCOMLIB_HPP)
#define LOOS_BCOMLIB_HPP


#include <loos.hpp>

/*
 * = Developer's Note =
 *
 * Some of the following is vestigial from the develpment of the BCOM method.
 * These interfaces and classes may go away or change substantially in
 * the next release or two, so caveat programmer...comments may also
 * be quite sparse...
 */



namespace Convergence {

  // Treats an AtomicGroup as a column vector and subtracts it from
  // each column of the matrix M.
  template<typename T>
  void subtractStructure(T& M, const loos::AtomicGroup& model) {
    std::vector<float> avg(model.size() * 3);
    uint k = 0;
    for (uint i=0; i<model.size(); ++i) {
      loos::GCoord c = model[i]->coords();
      avg[k++] = c.x();
      avg[k++] = c.y();
      avg[k++] = c.z();
    }
    
    for (uint i=0; i<M.cols(); ++i)
      for (uint j=0; j<M.rows(); ++j)
        M(j, i) -= avg[j];
  }

  // Computes the cosine content for a col-vector
  template<typename T>
  double cosineContent(T& V, const uint col) {
    double sum1 = 0;
    double sum2 = 0;

    uint m = V.rows();
    double k = (col+1) * M_PI / m;
    for (uint j=0; j<m; ++j) {
      sum1 += cos(k * j) * V(j, col);
      sum2 += V(j, col) * V(j, col);
    }
    double c = 2.0 * sum1 * sum1 / (sum2 * m);
    return(c);
  }
  

  
  /*
   * Various policies that determine how blocks are extracted and
   * averaged/aligned.  The idea is that they are really functors
   * which, given a vector<AtomicGroup> ensemble, will extract a
   * RealMatrix of coordinates where each AtomicGroup is a column
   * vector.  The appropriate processing (i.e. average subtraction) is
   * also performed by the functor.
   *
   * local_average, when set, means that the average of the ensemble
   * is used rather than the average passed to the constructor
   * (presumably, the average of the -entire- trajectory)
   */


  // Align to the passed structure
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


  // Do no aligning
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





  // Compute the PCA of an ensemble using the specified coordinate
  // extraction policy...
  //

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



  // Get just the RSVs (this is for cosine-content calculations)
  // given an extraction policy...
  //

  template<class ExtractPolicy>
  loos::RealMatrix rsv(std::vector<loos::AtomicGroup>& ensemble, ExtractPolicy& extractor) {

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

    // Correctly scale the eigenvalues
    for (uint j=0; j<W.rows(); ++j)
      W[j] = W[j] < 0 ? 0.0 : sqrt(W[j]);

    // Multiply eigenvectors by inverse eigenvalues
    for (uint i=0; i<C.cols(); ++i) {
      double konst = (W[i] > 0.0) ? (1.0/W[i]) : 0.0;

      for (uint j=0; j<C.rows(); ++j)
        C(j, i) *= konst;
    }

    W.reset();
    loos::RealMatrix Vt = loos::Math::MMMultiply(C, M, true, false);
    loos::RealMatrix V = loos::Math::transpose(Vt);
    return(V);
  }


}

#endif


// @endcond PACKAGES_INTERNAL
