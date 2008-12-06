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
using namespace loos;

#if defined(__linux__)
extern "C" {
  void sgesvd_(char*, char*, int*, int*, float*, int*, float*, float*, int*, float*, int*, float*, int*, int*);
}
#endif

typedef float svdreal;

typedef Math::Matrix<svdreal, Math::ColMajor> Matrix;

#define SVDFUNC  sgesvd_

struct Globals {
  string model_name, traj_name;
  string alignment_string, svd_string;
  greal alignment_tol;
  bool include_source;
  int terms;
  string prefix;
  string avg_name;
  uint dcdmin, dcdmax;
  string header;
};




Globals globals;

string findBaseName(const string& s) {
  string result;

  int n = s.find('.');
  result = (n <= 0) ? s : s.substr(0, n);

  return(result);
}


// svd [opts] pdb dcd

void parseOptions(int argc, char *argv[]) {

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("align,a", po::value<string>(&globals.alignment_string)->default_value("name == 'CA'"), "Selection to align with")
      ("svd,s", po::value<string>(&globals.svd_string)->default_value("!(segid == 'BULK' || segid == 'SOLV' || hydrogen)"), "Selection to calculate the SVD of")
      ("tolerance,t", po::value<greal>(&globals.alignment_tol)->default_value(1e-6), "Tolerance for iterative alignment")
      ("terms,T", po::value<int>(&globals.terms), "# of terms of the SVD to output")
      ("average,A", po::value<string>(&globals.avg_name), "Write out the average structure to this filename")
      ("prefix,p", po::value<string>(), "Prefix SVD output filenames with this string")
      ("range,r", po::value<string>(), "Range of frames from the trajectory to operate over")
      ("source,s", po::value<bool>(&globals.include_source)->default_value(false),"Write out the source conformation matrix");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&globals.model_name), "Model filename")
      ("traj", po::value<string>(&globals.traj_name), "Trajectory filename");

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
      cerr << generic;
      exit(-1);
    }

    if (vm.count("prefix"))
      globals.prefix = vm["prefix"].as<string>();
    else
      globals.prefix = findBaseName(globals.traj_name);

    if (vm.count("range")) {
      string rangespec = vm["range"].as<string>();
      int i = sscanf(rangespec.c_str(), "%u:%u", &globals.dcdmin, &globals.dcdmax);
      if (i != 2) {
        cerr << "Error - bad range given\n";
        exit(-1);
      }
    }
  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}


vector<XForm> doAlign(const AtomicGroup& subset, pTraj traj) {

  boost::tuple<vector<XForm>, greal, int> res = iterativeAlignment(subset, traj, globals.alignment_tol, 100);
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
  avgpdb.remarks().add(globals.header);
  ofstream ofs(globals.avg_name.c_str());
  if (!ofs)
    throw(runtime_error("Cannot open " + globals.avg_name));
  ofs << avgpdb;
}


// Calculates the transformed avg structure, then extracts the
// transformed coords from the DCD with the avg subtraced out...

Matrix extractCoords(const AtomicGroup& subset, const vector<XForm>& xforms, pTraj traj) {
  AtomicGroup avg = averageStructure(subset, xforms, traj);

  // Hook to get the avg structure if requested...
  if (globals.avg_name != "")
    writeAverage(avg);

  uint natoms = subset.size();
  AtomicGroup frame = subset.copy();
  uint n = globals.dcdmax - globals.dcdmin;
  uint m = natoms * 3;
  Matrix M(m, n);

  for (uint i=0; i<n; ++i) {
    traj->readFrame(i + globals.dcdmin);
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
  string header = invocationHeader(argc, argv);
  globals.header = header;

  parseOptions(argc, argv);

  // Need to address this...
  AtomicGroup model = createSystem(globals.model_name);
  pTraj ptraj = createTrajectory(globals.traj_name, model);
  
  // Fix max-range for DCD
  if (globals.dcdmax == 0)
    globals.dcdmax = ptraj->nframes();
  if (globals.dcdmin > ptraj->nframes() || globals.dcdmax > ptraj->nframes()) {
    cerr << "Invalid Trajectory range requested.\n";
    exit(-1);
  }

  AtomicGroup alignsub = selectAtoms(model, globals.alignment_string);
  AtomicGroup svdsub = selectAtoms(model, globals.svd_string);

  write_map(globals.prefix + ".map", svdsub);

  cerr << argv[0] << ": Aligning...\n";
  vector<XForm> xforms = doAlign(alignsub, ptraj);
  cerr << argv[0] << ": Extracting aligned coordinates...\n";
  Matrix A = extractCoords(svdsub, xforms, ptraj);
  f77int m = A.rows();
  f77int n = A.cols();
  f77int sn = m<n ? m : n;


  if (globals.include_source)
    writeAsciiMatrix(globals.prefix + "_A.asc", A, header);

  double estimate = m*m*sizeof(svdreal) + n*n*sizeof(svdreal) + m*n*sizeof(svdreal) + sn*sizeof(svdreal);
  cerr << argv[0] << ": Allocating space... (" << m << "," << n << ") for " << estimate/megabytes << "Mb\n";
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
  SVDFUNC(&jobu, &jobvt, &m, &n, A.get(), &lda, S.get(), U.get(), &ldu, Vt.get(), &ldvt, work, &lwork, &info);

  if (info > 0) {
    cerr << "Convergence error in dgesvd\n";
    exit(-3);
  } else if (info < 0) {
    cerr << "Error in " << info << "th argument to dgesvd\n";
    exit(-4);
  }
  cerr << argv[0] << ": Done!\n";

  writeAsciiMatrix(globals.prefix + "_U.asc", U, header);
  writeAsciiMatrix(globals.prefix + "_s.asc", S, header);
  writeAsciiMatrix(globals.prefix + "_V.asc", Vt, header, true);

  delete[] work;
}
