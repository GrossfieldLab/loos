/*
  big-svd

  Compute the SVD (PCA) of a large system/long trajectory.  This tool
  will use less memory than the svd tool since it does not compute all
  of the [unnecessary] right singular vectors.
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
namespace opts = loos::OptionsFramework;


const double KB = 1024.0;
const double MB = 1024 * KB;
const double GB = 1024 * MB;

// Don't include these classes in Doxygen
// @cond TOOLS_INTERNAL

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


class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() : write_source_matrix(false) { }

  void addGeneric(opts::po::options_description& o) {
    o.add_options()
      ("source", opts::po::value<bool>(&write_source_matrix)->default_value(write_source_matrix), "Write out source matrix");
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("source=%d") % write_source_matrix;
    return(oss.str());
  }

  bool write_source_matrix;
};
// @endcond







RealMatrix extractCoordinates(pTraj& traj, AtomicGroup& grp, const vector<uint>& indices) {
  uint m = grp.size() * 3;
  uint n = indices.size();

  RealMatrix A(m, n);
  vector<double> avg(m, 0.0);

  for (uint i=0; i<n; ++i) {
    traj->readFrame(indices[i]);
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




void writeMap(const string& fname, const AtomicGroup& grp) {
  ofstream fout(fname.c_str());

  if (!fout) {
    cerr << "Unable to open " << fname << " for output.\n";
    exit(-10);
  }

  pAtom pa;
  AtomicGroup::Iterator iter(grp);
  int i = 0;
  while (pa = iter())
    fout << i++ << "\t" << pa->id() << "\t" << pa->resid() << endl;

}


void normalizeRows(RealMatrix& A) {
  for (uint j=0; j<A.rows(); ++j) {
    double sum = 0.0;
    for (uint i=0; i<A.cols(); ++i)
      sum += A(j, i) * A(j, i);

    sum = sqrt(sum);
    for (uint i=0; i<A.cols(); ++i)
      A(j, i) /= sum;
  }
}


int main(int argc, char *argv[]) {

  string hdr = invocationHeader(argc, argv);
  
  opts::BasicOptions* bopts = new opts::BasicOptions;
  opts::BasicSelection* sopts = new opts::BasicSelection("name == 'CA'");
  opts::OutputPrefix* popts = new opts::OutputPrefix;
  opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(popts).add(tropts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  TrackStorage store;

  AtomicGroup model = tropts->model;
  AtomicGroup subset = selectAtoms(model, sopts->selection);
  pTraj traj = tropts->trajectory;

  string prefix = popts->prefix;
  vector<uint> indices = tropts->frameList();

  writeMap(prefix + ".map", subset);

  // Build AA'

  RealMatrix A = extractCoordinates(traj, subset, indices);
  cerr << boost::format("Coordinate matrix is %d x %d\n") % A.rows() % A.cols();
  store.allocate(A.rows() * A.cols());
  if (topts->write_source_matrix)
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

  // Multiply eigenvectors by inverse eigenvalues
  for (uint i=0; i<C.cols(); ++i) {
    double konst = (W[i] > 0.0) ? (1.0/W[i]) : 0.0;

    for (uint j=0; j<C.rows(); ++j)
      C(j, i) *= konst;
  }

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
