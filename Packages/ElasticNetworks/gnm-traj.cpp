/*
  gnm-traj

  Calculates a time-series of the first eigenvalue from a GNM calculated for each
  frame of a trajectory.

  See,
    Hall, B. A., Kaye, S. L., Pang, A., Perera, R. & Biggin, P. C. Characterization of protein conformational states by normal-mode frequencies. J Am Chem Soc 129, 11394–11401 (2007).

*/

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2013 Tod D. Romo
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


// This is the Kirchoff normalization constant (see Bahar, Atilgan,
// and Erman.  Folding & Design 2:173)
double normalization = 1.0;


DoubleMatrix kirchoff(AtomicGroup& group, const double cutoff) {
  int n = group.size();
  DoubleMatrix M(n, n);
  double r2 = cutoff * cutoff;


  for (int j=1; j<n; j++)
    for (int i=0; i<j; i++)
      if (group[i]->coords().distance2(group[j]->coords()) <= r2)
        M(i, j) = M(j, i) = -normalization;

  for (int j=0; j<n; j++) {
    double sum = 0;
    for (int i=0; i<n; i++) {
      if (i == j)
        continue;
      sum += M(j, i);
    }
    M(j, j) = -sum;
  }

  return(M);
}



int main(int argc, char *argv[]) 
{
  string hdr = invocationHeader(argc, argv);
  int k = 1;
  
  if (argc != 6) {
    cerr << "Usage- gnm-traj output-prefix model traj cutoff selection\n";
    exit(-1);
  }
  
  string prefix(argv[k++]);
  AtomicGroup model = createSystem(argv[k++]);
  pTraj traj = createTrajectory(argv[k++], model);
  double cutoff = strtod(argv[k++], 0);
  AtomicGroup subset = selectAtoms(model, argv[k++]);
  
  uint t = 0;

  

  uint n = subset.size();
  DoubleMatrix svals(traj->nframes(), 1);
  DoubleMatrix vecs(n, traj->nframes());

  cerr << "Progress- ";
  
  while (traj->readFrame()) {
    if (t % 100 == 0)
      cerr << '.';
    
    traj->updateGroupCoords(subset);
    DoubleMatrix K = kirchoff(subset, cutoff);
    boost::tuple<DoubleMatrix, DoubleMatrix, DoubleMatrix> result = svd(K);

    DoubleMatrix U = boost::get<0>(result);
    DoubleMatrix S = boost::get<1>(result);

    svals[t] = S[n - 2];
    for (uint i=0; i<n; ++i)
      vecs(i, t) = U(i, n-2);
    
    ++t;
  }

  cerr << endl;
  
  writeAsciiMatrix(prefix + "_s.asc", svals, hdr);
  writeAsciiMatrix(prefix + "_U.asc", vecs, hdr);

}


