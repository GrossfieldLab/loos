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
#include <boost/program_options.hpp>

using namespace std;
using namespace loos;
namespace po = boost::program_options;


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




// --------------------------- GLOBALS

vector<uint> indices;
string traj_name;
string model_name;
string prefix;
string selection;
bool write_source_matrix;



void parseArgs(int argc, char *argv[]) {
  
  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("range,r", po::value< vector<string> >(), "Range of frames from the trajectory to operate over")
      ("svd,s", po::value<string>(&selection)->default_value("name == 'CA'"), "Selection to calculate the SVD of")
      ("source,S", po::value<bool>(&write_source_matrix)->default_value(false), "Write out the source data matrix")
      ("prefix,p", po::value<string>(&prefix), "Output prefix");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value<string>(&traj_name), "Traj filename");


    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("traj", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("traj"))) {
      cout << "Usage- " << argv[0] << " [options] model traj output-prefix\n";
      cout << generic;
      exit(0);
    }

    if (vm.count("range")) {
      vector<string> ranges = vm["range"].as< vector<string> >();
      indices = parseRangeList<uint>(ranges);
    }

    if (vm.count("prefix"))
      prefix = vm["prefix"].as<string>();
    else
      prefix = findBaseName(traj_name);

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }

}




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




int main(int argc, char *argv[]) {

  string hdr = invocationHeader(argc, argv);
  parseArgs(argc, argv);

  TrackStorage store;

  AtomicGroup model = createSystem(model_name);
  AtomicGroup subset = selectAtoms(model, selection);
  pTraj traj = createTrajectory(traj_name, model);

  writeMap(prefix + ".map", subset);

  if (indices.empty())
    for (uint i=0; i<traj->nframes(); ++i)
      indices.push_back(i);

  // Build AA'

  RealMatrix A = extractCoordinates(traj, subset, indices);
  cerr << boost::format("Coordinate matrix is %d x %d\n") % A.rows() % A.cols();
  store.allocate(A.rows() * A.cols());
  if (write_source_matrix)
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
