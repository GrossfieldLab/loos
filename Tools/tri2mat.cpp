/*
  tri2mat
  
  Converts a triangular matrix into a full (square) matrix...
*/



/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo
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



#include <loos.hpp>

int main(int argc, char *argv[]) {

  string header = invocationHeader(argc, argv);
  loos::Matrix<float, Triangular> M;
  loos::readAsciiMatrix(argv[1], M);
  if (M.rows() == 0) {
    cerr << "Error- no rows in matrix.\n";
    exit(-1);
  }

  loos::writeTriangularAsFull(argv[2], M, header);
}
