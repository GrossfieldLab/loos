/*
  svd.cpp

  Computes the SVD for a trajectory.  Writes out the SVD as an
  OCTAVE-formatted text file.
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

#include <boost/program_options.hpp>

namespace po = boost::program_options;
using namespace std;
using namespace loos;

#if defined(__linux__)
extern "C" {
  void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
  void sgesvd_(char*, char*, int*, int*, float*, int*, float*, float*, int*, float*, int*, float*, int*, int*);
}
#endif


typedef double svdreal;
#define SVDFUNC  dgesvd_


typedef Math::Matrix<svdreal, Math::ColMajor> Matrix;


string model_name, traj_name;
string alignment_string, svd_string;
bool noalign;
greal alignment_tol;
bool include_source;
int terms;
string prefix;
string header;
vector<uint> indices;





void parseOptions(int argc, char *argv[]) {

  try {
    po::options_description generic("Allowed options", 120);
    generic.add_options()
      ("help,h", "Produce this help message")
      ("align,a", po::value<string>(&alignment_string)->default_value("name == 'CA'"), "Selection to align with")
      ("svd,s", po::value<string>(&svd_string)->default_value("!(segid == 'BULK' || segid == 'SOLV' || hydrogen)"), "Selection to calculate the SVD of")
      ("tolerance,t", po::value<greal>(&alignment_tol)->default_value(1e-6), "Tolerance for iterative alignment")
      ("noalign,n", po::value<bool>(&noalign)->default_value(false), "Do NOT align the frames of trajectory")
      ("terms,T", po::value<int>(), "# of terms of the SVD to output")
      ("prefix,p", po::value<string>(), "Prefix SVD output filenames with this string")
      ("range,r", po::value< vector<string> >(), "Range of frames from the trajectory to operate over")
      ("source,S", po::value<bool>(&include_source)->default_value(false),"Write out the source conformation matrix");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value<string>(&traj_name), "Trajectory filename");

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
      cerr << "Usage- svd [options] model-name trajectory-name\n";
      cerr << setprecision(2) << generic;
      exit(-1);
    }

    if (vm.count("terms"))
      terms = vm["terms"].as<int>();
    else
      terms = 0;

    if (vm.count("prefix"))
      prefix = vm["prefix"].as<string>();
    else
      prefix = findBaseName(traj_name);

    if (vm.count("range")) {
      vector<string> ranges = vm["range"].as< vector<string> >();
      indices = parseRangeList<uint>(ranges);
    }

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}


vector<XForm> doAlign(const AtomicGroup& subset, pTraj traj) {

  boost::tuple<vector<XForm>, greal, int> res = iterativeAlignment(subset, traj, indices, alignment_tol, 100);
  vector<XForm> xforms = boost::get<0>(res);
  greal rmsd = boost::get<1>(res);
  int iters = boost::get<2>(res);

  cerr << "Subset alignment with " << subset.size()
       << " atoms converged to " << rmsd << " rmsd after "
       << iters << " iterations.\n";

  return(xforms);
}


void writeAverage(const AtomicGroup& avg) {
  PDB avgpdb = PDB::fromAtomicGroup(avg);
  avgpdb.pruneBonds();
  avgpdb.remarks().add(header);

  string fname(prefix + "_avg.pdb");
  ofstream ofs(fname.c_str());
  if (!ofs)
    throw(runtime_error("Cannot open " + fname));
  ofs << avgpdb;
}


// Calculates the transformed avg structure, then extracts the
// transformed coords from the DCD with the avg subtraced out...

Matrix extractCoords(const AtomicGroup& subset, const vector<XForm>& xforms, pTraj traj) {

  AtomicGroup avg = averageStructure(subset, xforms, traj, indices);
  writeAverage(avg);

  uint natoms = subset.size();
  AtomicGroup frame = subset.copy();
  uint n = indices.size();
  uint m = natoms * 3;
  Matrix M(m, n);

  for (uint i=0; i<n; ++i) {
    traj->readFrame(indices[i]);
    traj->updateGroupCoords(frame);
    frame.applyTransform(xforms[i]);

    for (uint j=0; j<natoms; j++) {
      GCoord c = frame[j]->coords() - avg[j]->coords();
      M(j*3,i) = c.x();
      M(j*3+1,i) = c.y();
      M(j*3+2,i) = c.z();
    }
  }

  return(M);
}



void write_map(const string& fname, const AtomicGroup& grp) {
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
  header = invocationHeader(argc, argv);

  parseOptions(argc, argv);

  // Need to address this...
  AtomicGroup model = createSystem(model_name);
  pTraj ptraj = createTrajectory(traj_name, model);
  AtomicGroup svdsub = selectAtoms(model, svd_string);

  if (indices.empty())
    for (uint i=0; i<ptraj->nframes(); ++i)
      indices.push_back(i);

  write_map(prefix + ".map", svdsub);

  vector<XForm> xforms;
  if (noalign) {
    // Make noop xforms to prevent doing any alignment...
    cerr << argv[0] << ": SKIPPING ALIGNMENT\n";
    for (uint i=0; i<indices.size(); ++i)
      xforms.push_back(XForm());
  } else {
    AtomicGroup alignsub = selectAtoms(model, alignment_string);
    cerr << argv[0] << ": Aligning...\n";
    xforms = doAlign(alignsub, ptraj);   // Honors indices
  }

  cerr << argv[0] << ": Extracting coordinates...\n";
  Matrix A = extractCoords(svdsub, xforms, ptraj);   // Honors indices
  f77int m = A.rows();
  f77int n = A.cols();
  f77int sn = m<n ? m : n;


  if (include_source)
    writeAsciiMatrix(prefix + "_A.asc", A, header);

  double estimate = m*m*sizeof(svdreal) + n*n*sizeof(svdreal) + m*n*sizeof(svdreal) + sn*sizeof(svdreal);
  cerr << boost::format("%s: Allocating estimated %.2f MB for %d x %d SVD\n")
    % argv[0]
    % (estimate / megabytes)
    % m
    % n;

  char jobu = 'A', jobvt = 'A';
  f77int lda = m, ldu = m, ldvt = n, lwork= -1, info;
  svdreal prework[10], *work;

  Matrix U(m,m);
  Matrix S(sn,1);
  Matrix Vt(n,n);
  
  // First, request the optimal size of the work array...
  SVDFUNC(&jobu, &jobvt, &m, &n, A.get(), &lda, S.get(), U.get(), &ldu, Vt.get(), &ldvt, prework, &lwork, &info);
  if (info != 0) {
    cerr << "Error code from size request to dgesvd was " << info << endl;
    exit(-2);
  }

  lwork = (f77int)prework[0];
  estimate += lwork * sizeof(svdreal);
  cerr << argv[0] << ": SVD requests " << lwork << " extra space for a grand total of " << estimate / megabytes << "Mb\n";
  work = new svdreal[lwork];

  cerr << argv[0] << ": Calculating SVD...\n";
  Timer<WallTimer> timer;
  timer.start();
  SVDFUNC(&jobu, &jobvt, &m, &n, A.get(), &lda, S.get(), U.get(), &ldu, Vt.get(), &ldvt, work, &lwork, &info);
  timer.stop();
  cerr << argv[0] << ": Done!  Calculation took " << timeAsString(timer.elapsed()) << endl;

  if (info > 0) {
    cerr << "Convergence error in dgesvd\n";
    exit(-3);
  } else if (info < 0) {
    cerr << "Error in " << info << "th argument to dgesvd\n";
    exit(-4);
  }


  Math::Range orig(0,0);
  Math::Range Usize(m,m);
  Math::Range Ssize(sn,1);
  Math::Range Vsize(sn,n);

  if (terms > 0) {
    if (terms > m || terms > sn || terms > n) {
      cerr << "ERROR- The number of terms requested exceeds matrix dimensions.\n";
      exit(-1);
    }
    Usize = Math::Range(m, terms);
    Ssize = Math::Range(terms, 1);
    Vsize = Math::Range(terms, n);
  }

  cerr << argv[0] << ": Writing results...\n";
  writeAsciiMatrix(prefix + "_U.asc", U, header, orig, Usize);
  writeAsciiMatrix(prefix + "_s.asc", S, header, orig, Ssize);
  writeAsciiMatrix(prefix + "_V.asc", Vt, header, orig, Vsize, true);
  cerr << argv[0] << ": done!\n";

  delete[] work;
}
