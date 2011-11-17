/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2010, Tod D. Romo, Alan Grossfield
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

/** \addtogroup ENM
 *@{
 */

#if !defined(LOOS_ANM_LIB_HPP)
#define LOOS_ANM_LIB_HPP


#include "enm-lib.hpp"

namespace ENM {

  //! Anisotropic network model
  class ANM : public ElasticNetworkModel {
  public:
    ANM(SuperBlock* b) : ElasticNetworkModel(b) { prefix_ = "anm"; }

    void solve() {

      if (verbosity_ > 1)
        std::cerr << "Building hessian...\n";
      buildHessian();
      if (debugging_)
        loos::writeAsciiMatrix(prefix_ + "_H.asc", hessian_, meta_, false);

      loos::Timer<> t;
      if (verbosity_ > 0)
        std::cerr << "Computing SVD of hessian...\n";
      t.start();

      boost::tuple<loos::DoubleMatrix, loos::DoubleMatrix, loos::DoubleMatrix> result = svd(hessian_);
    
      t.stop();
      if (verbosity_ > 0)
        std::cerr << "SVD took " << loos::timeAsString(t.elapsed()) << std::endl;


      eigenvecs_ = boost::get<0>(result);
      eigenvals_ = boost::get<1>(result);
      rsv_ = boost::get<2>(result);

      loos::Math::reverseRows(eigenvals_);
      loos::Math::reverseColumns(eigenvecs_);
      loos::Math::reverseRows(rsv_);
    }


    //! Return the inverted hessian matrix
    loos::DoubleMatrix inverseHessian() {

      if (rsv_.rows() == 0)
        throw(std::logic_error("ANM::inverseHessian() called before ANM::solve()"));

      uint n = eigenvals_.rows();
      for (uint i=6; i<n; ++i) {
        double s = 1.0 / eigenvals_[i];
        for (uint j=0; j<n; ++j)
          rsv_(i, j) *= s;
      }

      loos::DoubleMatrix Hi = loos::Math::MMMultiply(rsv_, eigenvecs_, true, true);
      return(Hi);
    }


  private:
    loos::DoubleMatrix rsv_;

  };

};

#endif


/** @} */
