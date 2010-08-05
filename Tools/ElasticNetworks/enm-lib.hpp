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



#if defined(__linux__)
extern "C" {
  void dsygvx_(int*, char*, char*, char*, int*, double*, int*, double*, int*, double*, double*, int*, int*, double*, int*, double*, double*, int*, double*, int*, int*, int*, int*);
  void dpotrf_(char*, int*, double*, int*, int*);
  void dtrmm_(char*, char*, char*, char*, int*, int*, double*, double*, int*, double*, int*);
}
#endif



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



#endif
