/*
  MatrixIO.hpp
  Classes for reading and writing matrices in various formats...
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
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

#include <stdexcept>


#include <loos.hpp>
#include <MatrixIO.hpp>




void MatrixWriter::basic_write(const Wrapper& w, const string& tag, const uint m, const uint n, const bool trans, const uint maxcol, const uint maxrow) {
  uint i, j, k, s = m*n;
  uint nn = (maxcol == 0 || maxcol > n) ? n : maxcol;
  uint mm = (maxrow == 0 || maxrow > n) ? m : maxrow;
  ofstream *ofsp = 0;
  
  // Determine if we need to open a file or not...
  ostream *po = ofs;
  if (!po) {
    string fname = constructFilename(tag);
    ofsp = new ofstream(fname.c_str());
    if (!ofsp)
      throw(runtime_error("Unable to open file " + fname));
    po = ofsp;
  }

  OutputPreamble(po, tag, m, n, trans);
  for (j=0; j<mm; j++) {
    for (i=0; i<nn; i++) {
      double d;

      if (trans)
	k = j*n+i;
      else
	k = i*m+j;
      assert(k < s);
      d = w(k);
      OutputDatum(po, d);
    }
    OutputEOL(po);
  }
  OutputCoda(po);
  
  delete ofsp;
}
