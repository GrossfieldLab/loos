/*
  vsa

  (c) 2009 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry


  Usage:
    vsa subset environment radius model output_prefix

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



#include <loos.hpp>

#include <boost/format.hpp>
#include <boost/program_options.hpp>


using namespace std;
using namespace loos;
namespace po = boost::program_options;


typedef pair<uint,uint> Range;




const double normalization = 1.0;
double threshold = 1e-6;


DoubleMatrix hblock(const int i, const int j, const AtomicGroup& model, const double radius2) {

  DoubleMatrix B(3,3);
  GCoord u = model[i]->coords();
  GCoord v = model[j]->coords();
  GCoord d = v - u;

  double s = d.length2();
  if (s <= radius2) {

    for (int j=0; j<3; ++j)
      for (int i=0; i<3; ++i)
        B(i,j) = normalization * d[i]*d[j] / s;
  }

  return(B);
}



DoubleMatrix hessian(const AtomicGroup& model, const double radius) {
  
  int n = model.size();
  DoubleMatrix H(3*n,3*n);
  double r2 = radius * radius;

  for (int i=1; i<n; ++i) {
    for (int j=0; j<i; ++j) {
      DoubleMatrix B = hblock(i, j, model, r2);
      for (int x = 0; x<3; ++x)
        for (int y = 0; y<3; ++y) {
          H(i*3 + y, j*3 + x) = -B(y, x);
          H(j*3 + x, i*3 + y) = -B(x ,y);
        }
    }
  }

  // Now handle the diagonal...
  for (int i=0; i<n; ++i) {
    DoubleMatrix B(3,3);
    for (int j=0; j<n; ++j) {
      if (j == i)
        continue;
      
      for (int x=0; x<3; ++x)
        for (int y=0; y<3; ++y)
          B(y,x) += H(j*3 + y, i*3 + x);
    }

    for (int x=0; x<3; ++x)
      for (int y=0; y<3; ++y)
        H(i*3 + y, i*3 + x) = -B(y,x);
  }

  return(H);
}



DoubleMatrix submatrix(const DoubleMatrix& M, const Range& rows, const Range& cols) {
  uint m = rows.second - rows.first + 1;
  uint n = cols.second - cols.first + 1;

  DoubleMatrix A(m,n);
  for (uint i=0; i < n; ++i)
    for (uint j=0; j < m; ++j)
      A(j,i) = M(j+rows.first, i+cols.first);

  return(A);
}


// Sorts a vector of indices based on the first column of the bound matrix...
struct SortPredicate {
  SortPredicate(const DoubleMatrix& A) : M(A) { }
  bool operator()(const int i, const int j) { return(M[i]> M[j]); }

  const DoubleMatrix& M;
};


boost::tuple<DoubleMatrix, DoubleMatrix> eigenDecomp(DoubleMatrix& A, DoubleMatrix& B) {
  char jobvl = 'N';
  char jobvr = 'V';
  f77int fn = A.rows();
  f77int lda = fn;
  f77int ldb = fn;
  DoubleMatrix AR(fn,1);
  DoubleMatrix AI(fn,1);
  DoubleMatrix BETA(fn,1);
  double vl;
  f77int ldvl = 1;
  DoubleMatrix VR(fn,fn);
  f77int ldvr = fn;
  f77int lwork = 64 * fn;
  DoubleMatrix WORK(lwork, 1);

  f77int info;

  cout << "Calling dggev...\n";
  dggev_(&jobvl, &jobvr, &fn, A.get(), &lda, B.get(), &ldb, AR.get(), AI.get(), BETA.get(),
         &vl, &ldvl, VR.get(), &ldvr, WORK.get(), &lwork, &info);
  cout << "Result = " << info << endl;

  // Complex eigenvalues are set to 0
  DoubleMatrix eigvals(fn, 1);
  for (int i=0; i<fn; ++i) {
    if (AI[i] == 0.0)
      eigvals[i] = AR[i] / BETA[i];
    else
      eigvals[i] = 0.0;
  }

  // normalize
  for (int i=0; i<fn; ++i) {
    double norm = 0.0;
    for (int j=0; j<fn; ++j)
      norm += (VR(j,i) * VR(j,i));
    norm = sqrt(norm);
    for (int j=0; j<fn; ++j)
      VR(j,i) /= norm;
  }

  // Sort by eigenvalue
  vector<int> indices;
  for (int i=0; i<fn; ++i)
    indices.push_back(i);
  SortPredicate sp(eigvals);
  sort(indices.begin(), indices.end(), sp);

  // Now copy over in order...
  DoubleMatrix SD(fn, 1);
  DoubleMatrix SU(fn, fn);
  for (int i=0; i<fn; ++i) {
    SD(i,0) = eigvals[indices[i]];
    for (int j=0; j<fn; ++j)
      SU(j,i) = VR(j,indices[i]);
  }

  boost::tuple<DoubleMatrix, DoubleMatrix> result(SD, SU);
  return(result);
}


void assignMasses(AtomicGroup& grp, const string& name) {
  ifstream ifs(name.c_str());
  if (!ifs) {
    cerr << "Warning- no masses will be used\n";
    return;
  }

  double mass;
  string pattern;
  bool no_masses_set = true;

  while (ifs >> pattern >> mass) {
    string selection("name =~ '");
    selection += pattern + "'";
    Parser parser(selection);
    KernelSelector selector(parser.kernel());
    AtomicGroup subset = grp.select(selector);
    if (subset.empty())
      continue;
    no_masses_set = false;
    cerr << "Assigning " << subset.size() << " atoms of type " << pattern << " to mass " << mass << endl;
    for (AtomicGroup::iterator i = subset.begin(); i != subset.end(); ++i)
      (*i)->mass(mass);
  }

  if (no_masses_set)
    cerr << "WARNING- no masses were assigned!\n";
}


DoubleMatrix getMasses(const AtomicGroup& grp) {
  uint n = grp.size();

  DoubleMatrix M(3*n,3*n);
  for (uint i=0, k=0; i<n; ++i, k += 3) {
    M(k,k) = grp[i]->mass();
    M(k+1,k+1) = grp[i]->mass();
    M(k+2,k+2) = grp[i]->mass();
  }

  return(M);
}


void showSize(const string& s, const DoubleMatrix& M) {
  cout << s << M.rows() << " x " << M.cols() << endl;
}


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  uint k=1;

  string subsel(argv[k++]);
  string envsel(argv[k++]);
  double rad = strtod(argv[k++], 0);
  
  AtomicGroup model = createSystem(argv[k++]);
  string prefix(argv[k++]);
  assignMasses(model, argv[k++]);

  AtomicGroup subset = selectAtoms(model, subsel);
  AtomicGroup environment = selectAtoms(model, envsel);
  AtomicGroup composite = subset + environment;

  cerr << "Subset size is " << subset.size() << endl;
  cerr << "Environment size is " << environment.size() << endl;

  DoubleMatrix H = hessian(composite, rad);
  showSize("H: ", H);
  // Now, burst out the subparts...
  uint l = subset.size() * 3;

  uint n = H.cols() - 1;
  DoubleMatrix Hss = submatrix(H, Range(0,l-1), Range(0,l-1));
  showSize("Hss: ", Hss);

  DoubleMatrix Hee = submatrix(H, Range(l, n), Range(l, n));
  showSize("Hee: ", Hee);

  DoubleMatrix Hse = submatrix(H, Range(0,l-1), Range(l, n));
  showSize("Hse: ", Hse);

  DoubleMatrix Hes = submatrix(H, Range(l, n), Range(0, l-1));
  showSize("Hes: ", Hes);

  cerr << "Inverting environment hessian...\n";
  DoubleMatrix Heei = Math::invert(Hee);
  cerr << "Computing Hssp...\n";
  DoubleMatrix Hssp = Hss - Hse * Heei * Hes;

  DoubleMatrix Ms = getMasses(subset);
  showSize("Ms: ", Ms);
  DoubleMatrix Me = getMasses(environment);
  showSize("Me: ", Me);

  cerr << "Computing Msp...\n";
  DoubleMatrix Msp = Ms + Hse * Heei * Me * Heei * Hes;

  cout << "Running eigendecomp of " << Hssp.rows() << " x " << Hssp.cols() << " matrix\n";
  boost::tuple<DoubleMatrix, DoubleMatrix> eigenpairs = eigenDecomp(Hssp, Msp);
  cout << "done\n";
  DoubleMatrix Ds = boost::get<0>(eigenpairs);
  DoubleMatrix Us = boost::get<1>(eigenpairs);

  writeAsciiMatrix(prefix + "_Ds.asc", Ds, hdr);
  writeAsciiMatrix(prefix + "_Us.asc", Us, hdr);

  

//   DoubleMatrix Ue = -Heei * Hes * Us;
//   X = Heei * Me * Heei * Hes * Us;
//   Y = -Ds;
//   DoubleMatrix De = mmult(Y, X, true, false);
//   writeAsciiMatrix(prefix + "_De.asc", De, hdr, true);
//   writeAsciiMatrix(prefix + "_Ue.asc", Ue, hdr);
}
