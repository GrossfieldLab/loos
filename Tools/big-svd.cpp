/*
  big-svd

  Compute the SVD (PCA) of a large system/long trajectory.  This tool
  should use less memory than the svd tool does.
*/



/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo
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
#include <boost/format.hpp>

using namespace std;
using namespace loos;


const double KB = 1024.0;
const double MB = 1024 * KB;
const double GB = 1024 * MB;


struct TrackStorage {
  TrackStorage() : storage(0) { }
  
  void allocate(ulong n) {

    n *= sizeof(float);
    storage += n;
    cerr << boost::format("Allocated %s for a total of %s memory\n")
      % memory(n)
      % memory(storage);
  }

  void free(const ulong n) {
    storage -= n;
  }


  string memory(const ulong n) {
    double val = static_cast<double>(n);
    string units;

    if (val >= GB) {
      val /= GB;
      units = "GB";
    } else if (val >= MB) {
      val /= MB;
      units = "MB";
    } else if (val >= KB) {
      val /= KB;
      units = "KB";
    } else
      units = "Bytes";

    ostringstream oss;
    oss << boost::format("%.2f %s") % val % units;
    return(oss.str());
  }

  ulong storage;
};




RealMatrix extractCoordinates(pTraj& traj, AtomicGroup& grp) {
  uint m = grp.size() * 3;
  uint n = traj->nframes();

  RealMatrix A(m, n);
  vector<double> avg(m, 0.0);

  for (uint i=0; i<n; ++i) {
    traj->readFrame(i);
    traj->updateGroupCoords(grp);
    for (uint j=0; j<static_cast<uint>(grp.size()); ++j) {
      GCoord c = grp[j]->coords();
      A(3*j, i) = c.x();
      avg[3*j] += c.x();

      A(3*j+1, i) = c.y();
      avg[3*j+1] += c.y();

      A(3*j+2, i) = c.z();
      avg[3*j+2] += c.z();
    }
  }

  for (uint j=0; j<m; ++j)
    avg[j] /= n;

  for (uint i=0; i<n; ++i)
    for (uint j=0; j<m; ++j)
      A(j, i) -= avg[j];

  return(A);
}





int main(int argc, char *argv[]) {
  if (argc == 1) {
    cerr << "Usage- big-svd selection model traj prefix\n";
    exit(0);
  }

  TrackStorage store;

  string hdr = invocationHeader(argc, argv);
  int k = 1;
  string selection(argv[k++]);
  string modelname(argv[k++]);
  string trajname(argv[k++]);
  string prefix(argv[k++]);

  AtomicGroup model = createSystem(modelname);
  AtomicGroup subset = selectAtoms(model, selection);
  pTraj traj = createTrajectory(trajname, model);

  // Build AA'

  RealMatrix A = extractCoordinates(traj, subset);
  cerr << boost::format("Coordinate matrix is %d x %d\n") % A.rows() % A.cols();
  store.allocate(A.rows() * A.cols());
  writeAsciiMatrix(prefix + "_A.asc", A, hdr);


  store.allocate(A.rows() * A.rows());
  cerr << "Multiplying transpose...\n";
  RealMatrix C = MMMultiply(A, A, false, true);
  cerr << "Done!\n";

  // Compute [U,D] = eig(C)

  char jobz = 'V';
  char uplo = 'L';
  f77int n = A.rows();
  f77int lda = n;
  float dummy;
  RealMatrix W(n, 1);
  f77int lwork = -1;
  f77int info;

  cerr << "Calling ssyev to get work size...\n";

  ssyev_(&jobz, &uplo, &n, C.get(), &lda, W.get(), &dummy, &lwork, &info);
  if (info != 0) {
      cerr << boost::format("ssyev failed with info = %d\n") % info;
      exit(-10);
  }
   
  lwork = static_cast<f77int>(dummy);
  store.allocate(lwork);
  float *work = new float[lwork+1];

  cerr << "Calling ssyev for eigendecomp...\n";
  ssyev_(&jobz, &uplo, &n, C.get(), &lda, W.get(), work, &lwork, &info);
  if (info != 0) {
      cerr << boost::format("ssyev failed with info = %d\n") % info;
      exit(-10);
  }
  cerr << "Finished!\n";
  
  reverseColumns(C);
  writeAsciiMatrix(prefix + "_U.asc", C, hdr);

  // D = sqrt(D);  Scale eigenvectors...
  for (uint j=0; j<W.rows(); ++j)
    W[j] = W[j] < 0 ? 0.0 : sqrt(W[j]);

  reverseRows(W);
  writeAsciiMatrix(prefix + "_s.asc", W, hdr);
  W.reset();

  store.free(W.rows() * W.cols());

  store.allocate(A.cols() * A.rows());
  cerr << "Multiplying to get RSVs...\n";
  RealMatrix Vt = MMMultiply(C, A, true, false);
  cerr << "Done!\n";
  C.reset();
  A.reset();
  writeAsciiMatrix(prefix + "_V.asc", Vt, hdr, true);

}
