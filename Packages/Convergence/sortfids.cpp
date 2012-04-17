/*
  sortfids

  Structural histogram, a la Lyman &
  Zuckerman, Biophys J (2006) 91:164-172

  Usage- sortfids model selection fids hist newfidname

  Sorts fiducials based on decreasing histogram bin population

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



#include <loos.hpp>

using namespace std;
using namespace loos;


typedef vector<AtomicGroup>                   vGroup;
typedef Math::Matrix<double, Math::RowMajor>  Matrix;

// @cond TOOLS_INTERNAL



string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\tSorts fiducial structures based on histogram bin population\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tGiven a set of fiducials for a structural histogram, sort them based\n"
    "on bin population.  This can be useful when using fidpick, which selects\n"
    "fiducials based on distance rather than bin probability.\n"
    "\n"
    "SEE ALSO\n"
    "\tfidpick\n";

  return(msg);
}


// Sorts based on 3rd col of a matrix

struct Adapter {
  Adapter(const Matrix& M) : A(M) { }

  uint size() const { return(A.rows()); }
  const double& operator[](const uint i) const {
    return(A(i,2));
  }

  const Matrix& A;

};

// @endcond

int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  
  if (argc != 6) {
    cerr << "Usage- sortfids model sel fids hist newfids\n";
    cerr << fullHelpMessage();
    exit(-1);
  }

  int opti = 1;

  AtomicGroup model = createSystem(argv[opti++]);
  string selection(argv[opti++]);
  AtomicGroup subset = selectAtoms(model, selection);
  subset.renumber();

  pTraj fids = createTrajectory(argv[opti++], subset);
  vGroup fiducials;
  readTrajectory(fiducials, subset, fids);

  Matrix M;
  readAsciiMatrix(argv[opti++], M);
  vector<uint> indices = sortedIndex<Adapter, DescendingSort<Adapter> >(Adapter(M));

  vGroup sorted;
  
  Matrix A(M.rows(), M.cols());
  for (uint i=0; i<M.rows(); ++i) {
    sorted.push_back(fiducials[indices[i]]);
    A(i,0) = i;
    A(i,1) = M(indices[i], 1);
    A(i,2) = M(indices[i], 2);
  }

  DCDWriter(argv[opti++], sorted, hdr);
  writeAsciiMatrix(cout, A, hdr);
}
