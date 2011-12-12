/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2010 Tod D. Romo
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


#if !defined(LOOS_VSA_LIB_HPP)
#define LOOS_VSA_LIB_HPP

#include <loos.hpp>
#include "enm-lib.hpp"

#if defined(__linux__) || defined(__CYGWIN__) || defined(__FreeBSD__)
extern "C" {
  void dsygvx_(int*, char*, char*, char*, int*, double*, int*, double*, int*, double*, double*, int*, int*, double*, int*, double*, double*, int*, double*, int*, int*, int*, int*);
  void dpotrf_(char*, int*, double*, int*, int*);
  void dtrmm_(char*, char*, char*, char*, int*, int*, double*, double*, int*, double*, int*);
}
#endif


namespace ENM {


  //! Vibrational subsystem analysis ENM
  /**
   *  References:
   *    - <a href="http://dx.doi.org/10.1063/1.3013558">Woodcock et al, J Chem Phys (2008) 129:214109</a>
   *    - <a href="http://dx.doi.org/10.1063/1.3141022">Haffner & Zheng, J Chem Phys (2009) 130:194111</a>
   *
   * The VSA class expects that the list of nodes contained in the
   * passed SuperBlock instance represents the combined system,
   * i.e. subsystem and environment.  The first \a subn nodes are the
   * subsystem.
   */
  class VSA : public ElasticNetworkModel {
  public:

    //! Constructor for VSA without masses
    /**
     * Arguments:
     * \arg \c blocker Determines how the Hessian is built
     * \arg \c subn The number of nodes in the subsystem
     */
    VSA(SuperBlock* blocker, const uint subn) : 
      ElasticNetworkModel(blocker),
      subset_size_(subn)
    { prefix_ = "vsa"; }

    //! Constructor for VSA with masses
    /**
     * Arguments:
     * \arg c blocker Determines how the Hessian is built
     * \arg c subn The number of nodes in the subsystem
     * \arg c M Diagonal 3N x 3N matrix of node masses
     */
    VSA(SuperBlock* blocker, const uint subn, const loos::DoubleMatrix& M) :
      ElasticNetworkModel(blocker),
      subset_size_(subn),
      masses_(M)
    { prefix_ = "vsa"; }



    void solve();
  
    //! Sets the mass matrix and determines what kind of VSA calc to use
    /**
     * Setting the mass matrix to an initialized matrix implies that VSA
     * will use the mass-VSA version.  On the other hand, setting the
     * matrix to a default, uninitialized matrix will switch to the
     * mass-less VSA:
     \code
     vsa.setMasses(DoubleMatrix());
     \endcode
    */
    void setMasses(const loos::DoubleMatrix& M) {
      masses_ = M;
    };


    //! Free up internal storage...
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


};


#endif



/** @} */
