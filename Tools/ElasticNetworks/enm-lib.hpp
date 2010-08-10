/*
  enm-lib

  (c) 2010 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry


  Common code for the ENM suite

*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009 Tod D. Romo
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



#if !defined(ENMLIB_HPP)
#define ENMLIB_HPP

#include <loos.hpp>
#include "hessian.hpp"



typedef std::pair<uint,uint> Range;

loos::DoubleMatrix submatrix(const loos::DoubleMatrix& M, const Range& rows, const Range& cols);

void normalizeColumns(loos::DoubleMatrix& A);


// Map masses from one group onto another...  Minimal error checking...
void copyMasses(loos::AtomicGroup& target, const loos::AtomicGroup& source);



// Copy the masses from a PSF onto a group
void massFromPSF(loos::AtomicGroup& grp, const std::string& name);

// The masses are stored in the occupancy field of a PDB...
void massFromOccupancy(loos::AtomicGroup& grp);


// Build the 3n x 3n diagonal mass matrix for a group
loos::DoubleMatrix getMasses(const loos::AtomicGroup& grp);





class ElasticNetworkModel {
public:
  ElasticNetworkModel(SuperBlock* blocker) : blocker_(blocker), name_("ENM"), prefix_(""), meta_(""), debugging_(false), verbosity_(0) { }
  virtual ~ElasticNetworkModel() { }

  // Should we allow this?
  void setSuperBlockFunction(SuperBlock* p) { blocker_ = p; }

  virtual void solve() =0;



  void prefix(const std::string& s) { prefix_ = s; }
  std::string prefix() const { return(prefix_); }

  void meta(const std::string& s) { meta_ = s; }
  std::string meta() const { return(meta_); }

  void debugging(const bool b) { debugging_ = b; }
  bool debugging() const { return(debugging_); }

  void verbosity(const int i) { verbosity_ = i; }
  int verbosity() const { return(verbosity_); }

  SpringFunction::Params setParams(const SpringFunction::Params& v) {
    return(blocker_->setParams(v));
  }

  bool validParams() const { return(blocker_->validParams()); }

  const loos::DoubleMatrix& eigenvectors() const { return(eigenvecs_); }
  const loos::DoubleMatrix& eigenvalues() const { return(eigenvals_); }
  const loos::DoubleMatrix& hessian() const { return(hessian_); }

protected:
  

  void buildHessian() {
  
    uint n = blocker_->size();
    loos::DoubleMatrix H(3*n,3*n);

    for (uint i=1; i<n; ++i) {
      for (uint j=0; j<i; ++j) {
        loos::DoubleMatrix B = blocker_->block(i, j);
        for (uint x = 0; x<3; ++x)
          for (uint y = 0; y<3; ++y) {
            H(i*3 + y, j*3 + x) = -B(y, x);
            H(j*3 + x, i*3 + y) = -B(x ,y);
          }
      }
    }

    // Now handle the diagonal...
    for (uint i=0; i<n; ++i) {
      loos::DoubleMatrix B(3,3);
      for (uint j=0; j<n; ++j) {
        if (j == i)
          continue;
      
        for (uint x=0; x<3; ++x)
          for (uint y=0; y<3; ++y)
            B(y,x) += H(j*3 + y, i*3 + x);
      }

      for (uint x=0; x<3; ++x)
        for (uint y=0; y<3; ++y)
          H(i*3 + y, i*3 + x) = -B(y,x);
    }

    hessian_ = H;
  }


protected:
  // Arguably, some of the following should be private rather than
  // protected...  But for now, we'll just cheat and make 'em all
  // protected...   Nyah, nyah!
  SuperBlock* blocker_;
  std::string name_;
  std::string prefix_;
  std::string meta_;
  bool debugging_;
  int verbosity_;

  loos::DoubleMatrix eigenvecs_;
  loos::DoubleMatrix eigenvals_;

  loos::DoubleMatrix hessian_;
  
};


#endif
