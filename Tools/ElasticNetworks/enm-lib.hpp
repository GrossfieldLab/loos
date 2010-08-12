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


/** \defgroup ENM Elastic Network Models
 *@{
 */

#if !defined(ENMLIB_HPP)
#define ENMLIB_HPP

#include <loos.hpp>
#include "hessian.hpp"




// -------------------------------------
// Support routines & types

//! Specify a range of columns/rows [first,second)
typedef std::pair<uint,uint> Range;

//! Extracts a submatrix
loos::DoubleMatrix submatrix(const loos::DoubleMatrix& M, const Range& rows, const Range& cols);

//! Normalizes each column as a column-vector
void normalizeColumns(loos::DoubleMatrix& A);


//! Map masses from one group onto another...  Minimal error checking...
void copyMasses(loos::AtomicGroup& target, const loos::AtomicGroup& source);



//! Copy the masses from a PSF onto a group
void massFromPSF(loos::AtomicGroup& grp, const std::string& name);

//! The masses are stored in the occupancy field of a PDB...
void massFromOccupancy(loos::AtomicGroup& grp);


//! Build the 3n x 3n diagonal mass matrix for a group
loos::DoubleMatrix getMasses(const loos::AtomicGroup& grp);


// -------------------------------------



// This class describes the interface for all ENMs...
// To instantiate, you must pass a SuperBlock which determines how the
// hessian is built.

//! Interface for all ENMs
class ElasticNetworkModel {
public:
  ElasticNetworkModel(SuperBlock* blocker) : blocker_(blocker), name_("ENM"), prefix_(""), meta_(""), debugging_(false), verbosity_(0) { }
  virtual ~ElasticNetworkModel() { }

  // Should we allow this?
  void setSuperBlockFunction(SuperBlock* p) { blocker_ = p; }

  // This computes the hessian and solves for the eigenpairs
  virtual void solve() =0;

  // filename prefix when we have to write something out
  void prefix(const std::string& s) { prefix_ = s; }
  std::string prefix() const { return(prefix_); }

  // Any metadata that gets added to matrices written out
  void meta(const std::string& s) { meta_ = s; }
  std::string meta() const { return(meta_); }

  // Debugging flag (generally means write out all intermediate matrices)
  void debugging(const bool b) { debugging_ = b; }
  bool debugging() const { return(debugging_); }

  // How wordy are we?
  void verbosity(const int i) { verbosity_ = i; }
  int verbosity() const { return(verbosity_); }

  // -----------------------------------------------------
  // The following forward to the contained SuperBlock...
  SpringFunction::Params setParams(const SpringFunction::Params& v) {
    return(blocker_->setParams(v));
  }

  bool validParams() const { return(blocker_->validParams()); }

  uint paramSize() const { return(blocker_->paramSize()); }
  // -----------------------------------------------------


  // Accessors for eigenpairs and hessian
  const loos::DoubleMatrix& eigenvectors() const { return(eigenvecs_); }
  const loos::DoubleMatrix& eigenvalues() const { return(eigenvals_); }
  const loos::DoubleMatrix& hessian() const { return(hessian_); }



protected:
  

  // It is not expected that subclasses will want to override this...
  // Uses the contained SuperBlock to build a hessian
  void buildHessian();
  

protected:
  // Arguably, some of the following should be private rather than
  // protected...  But for now, we'll just cheat and make 'em all
  // protected..
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


/** @} */
