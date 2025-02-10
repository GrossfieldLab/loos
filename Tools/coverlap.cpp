/*
  coverlap

  (c) 2009 Tod D. Romo, Grossfield Lab
           Department of Biochemistry
           University of Rochster School of Medicine and Dentistry


  Covariance overlap between ENM and/or PCA results (i.e. eigenpairs)
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

string lefts_name, leftU_name, rights_name, rightU_name;
bool left_is_enm;
bool right_is_enm;
bool square_left;
bool square_right;
bool scale_power;
uint number_of_modes;
double lscale;
double rscale;
uint subspace_size;
uint ntries;
uint seed = 0;

uint skip;


void fullHelpMessage() {
  //string msg =
  cout <<  "\n"
  "SYNOPSIS\n"
  "\n"
  "Calculate the covariance overlap between two eigendecompositions\n"
  "\n"
  "DESCRIPTION\n"
  "\n"
  "Coverlap compares two sets of eigendompositions (eigenvalues and\n"
  "eigenvectors) using the metric described by Hess. 2002, Phys Rev E\n"
  "(10.1103/PhysRevE.65.031910).\n"
  "\n"
  "In simplest terms this is a cumulative dot product between two\n"
  "eigenvectors that is weighted relatively by each associated eigenvalue.\n"
  "The nomenclature of the input is broken down as:\n"
  "\tls - left side eigenvalue file (filename_s.asc if made with loos)\n"
  "\tlU - left side eigenvector file (filename_U.asc if made with loos)\n"
  "\trs - right side eigenvalue file (filename_s.asc if made with loos)\n"
  "\trU - right side eigenvector file (filename_U.asc if made with loos)\n"
  "\t\n"
  "\n"
  "The results commonly compared are simulation PCA's and/or ENM-NMA's.\n"
  "ENM eigenpairs are handled differently from PCA eigenpairs, you must specify\n"
  "which sides are ENM results.  This is done with the -e1 for ENM results on\n"
  "the left side and -E1 for ENM results on the right side. Additionally, PCA \n"
  "eigenpairs can be real eigenpairs or they can come from an SVD, in which \n"
  "case the 'eigenvalues' must be squared. This is done with -s1 to square the\n"
  "left side eigenvalues and -S1 to square the right side eigenvalues.\n"
  "Finally, when comparing ENM and PCA you will probably want to scale the\n"
  "eigenvalues such that the total power on each side are comparable.\n"
  "The -p1 option does this.\n"
  "\n"
  "In addition a non-weighted overlap between eigenvectors can be computed.\n"
  "Since the eigenvectors are complete basis sets this will always sum to 1\n"
  "Therefore, a subset of the modes should be used.  The number of modes to\n"
  "use with this calculation is specified with the -u## argument.\n"
  "\n"
  "Another option is to compute the Z-score of a covariance overlap result.\n"
  "This is done using a boot-strap approach to scramble eigenpairs and compute\n"
  "a new overlap.  this is done -z=### times, which is used to calculate the\n"
  "average and standard deviation for the Z-score.  For more details see:\n"
  "\n"
  "  Romo, T. D., and Grossfield, A. Block covariance overlap method and \n"
  "  convergence in molecular dynamics simulation, J. Chem. Theor. Comput.,\n"
  "  2011, 7, 2464-2472\n"
  "\n"
  "  Leioatts, N., Romo, T. D., and Grossfield, A., Elastic network models \n"
  "  are robust to variations in formalism, J. Chem. Theor. Comput., 2012, \n"
  "  8,  2424-2434\n"
  "\n"
  "\n"
  "EXAMPLES\n"
  "\n"
  "\tcoverlap -e1 -S1 -p1 -u50 anm_s.asc anm_U.asc pca_s.asc pca_U.asc\n"
  "This computes the covariance overlap between an ANM result (the left side)\n"
  "and a PCA (the right side) that came from an SVD (right side).  On the right\n"
  "the singular values are squared (to make them eigenvalues) and they are\n"
  "scaled to match the ANM eigenvalues (-p1).  Finally, a subspace overlap using\n"
  "the first 50 modes is also computed.\n"
  "\n"
  "\tcoverlap -e1 -p1 -u50 anm_s.asc anm_U.asc pca_s.asc pca_U.asc\n"
  "The same as the above example, but here the PCA came from an eigendecomp,\n"
  "so the eigenvalues used are real eigenvalues and do not need to be squared.\n"
  "\n"
  "\tcoverlap -e1 -E1 -u25 anm_s.asc anm_U.asc vsa_s.asc vsa_U.asc\n"
  "This computes the covariance overlap between an ANM and a VSA model.\n"
  "No scaling is applied to either side.  The subspace overlap using the\n"
  "first 25 modes is also computed.\n"
  "\n"
  "\tcoverlap -e1 -E1 -u25 -k 1.234 anm_s.asc anm_U.asc vsa_s.asc vsa_U.asc\n"
  "The same as the above example, but here 1.234 is used to scale the\n"
  "ANM eigenvalues.\n"
  "\n";
}



void parseArgs(int argc, char *argv[]) {

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("fullhelp", "Get extended help")
      ("skip,i", po::value<uint>(&skip)->default_value(6), "# of eigenvalues to skip for ENM")
      ("left_enm,e", po::value<bool>(&left_is_enm)->default_value(false), "Left side contains ENM results")
      ("right_enm,E", po::value<bool>(&right_is_enm)->default_value(false), "Right side contains ENM results")
      ("square_left,s", po::value<bool>(&square_left)->default_value(false), "Square left side (assumes PCA)")
      ("square_right,S", po::value<bool>(&square_right)->default_value(false), "Square right side (assumes PCA)")
      ("power,p", po::value<bool>(&scale_power)->default_value(false), "Scale the eigenvalue power of the right side to the left")
      ("modes,m", po::value<uint>(&number_of_modes)->default_value(0), "Number of modes to compare...  0 = all")
      ("left_scale,k", po::value<double>(&lscale)->default_value(1.0), "Scale left eigenvalues by this constant")
      ("right_scale,K", po::value<double>(&rscale)->default_value(1.0), "Scale right eigenvalues by this constant")
      ("subspace,u", po::value<uint>(&subspace_size)->default_value(25), "# of modes to use for the subspace overlap (0 = same as covariance)")
      ("zscore,z", po::value<uint>(&ntries)->default_value(0), "Use z-score (sets number of repeats)")
      ("seed", po::value<uint>(&seed)->default_value(0), "Seed for random number generator (0 = auto)");


    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("ls", po::value<string>(&lefts_name), "Left eigenvalues")
      ("lu", po::value<string>(&leftU_name), "Left eigenvector")
      ("rs", po::value<string>(&rights_name), "Right eigenvector")
      ("ru", po::value<string>(&rightU_name), "Right eigenvector");


    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("ls", 1);
    p.add("lu", 1);
    p.add("rs", 1);
    p.add("ru", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || vm.count("fullhelp") || !(vm.count("ls") && vm.count("lu") && vm.count("rs") && vm.count("ru"))) {
      cerr << "Usage- " << argv[0] << " [options] ls lU rs rU >output\n";
      cerr << generic;
      if (vm.count("fullhelp"))
        fullHelpMessage();
      exit(-1);
    }

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}

typedef boost::tuple<DoubleMatrix, DoubleMatrix>   DMDuple;


DMDuple transformENM(const DoubleMatrix& S, const DoubleMatrix& U, const uint nmodes) {
  DoubleMatrix SS(nmodes, 1);
  DoubleMatrix UU(U.rows(), nmodes);

  for (uint i=skip; i<nmodes+skip; ++i) {
    SS[i-skip] = 1.0 / S[i];
    for (uint j=0; j<U.rows(); ++j)
      UU(j,i-skip) = U(j,i);
  }

  return(DMDuple(SS,UU));
}


DMDuple firstColumns(const DoubleMatrix& S, const DoubleMatrix& U, const uint nmodes) {
  DoubleMatrix SS(nmodes, 1);
  DoubleMatrix UU(U.rows(), nmodes);

  for (uint i=0; i<nmodes; ++i) {
    SS[i] = i < S.rows() ? S[i] : 0.0;

    for (uint j=0; j<U.rows(); ++j)
      UU(j, i) = U(j, i);
  }

  return(DMDuple(SS,UU));
}



DoubleMatrix scalePower(const DoubleMatrix& A, const DoubleMatrix& B) {

  double sumB = 0.0;
  double sumA = 0.0;
  for (uint j=0; j<B.rows(); ++j) {
    sumB += B[j];
    sumA += A[j];
  }

  double scale = sumA / sumB;
  cerr << "Scale factor = " << scale << endl;
  DoubleMatrix E(B.rows(), 1);
  for (uint j=0; j<B.rows(); ++j)
    E[j] = B[j] * scale;

  return(E);
}




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  parseArgs(argc, argv);

  cerr << "Reading left side matrices...\n";
  DoubleMatrix lS;
  readAsciiMatrix(lefts_name, lS);
  DoubleMatrix lU;
  readAsciiMatrix(leftU_name, lU);
  cerr << boost::format("Read in %d x %d eigenvectors...\n") % lU.rows() % lU.cols();
  cerr << boost::format("Read in %d eigenvalues...\n") % lS.rows();

  cerr << "Reading in right side matrices...\n";
  DoubleMatrix rS;
  readAsciiMatrix(rights_name, rS);
  DoubleMatrix rU;
  readAsciiMatrix(rightU_name, rU);
  cerr << boost::format("Read in %d x %d eigenvectors...\n") % rU.rows() % rU.cols();
  cerr << boost::format("Read in %d eigenvalues...\n") % rS.rows();

  if (number_of_modes == 0) {
    number_of_modes = lS.rows() > rS.rows() ? lS.rows() : rS.rows();
    if (left_is_enm || right_is_enm)
      number_of_modes -= skip;
  }


  if (subspace_size > number_of_modes) {
    cerr << "ERROR- subspace size cannot exceed number of modes for covariance overlap\n";
    exit(-1);
  }

  DoubleMatrix lSS;
  DoubleMatrix lUU;
  if (left_is_enm) {
    DMDuple res = transformENM(lS, lU, number_of_modes);
    lSS = boost::get<0>(res);
    lUU = boost::get<1>(res);
  } else {
    DMDuple res = firstColumns(lS, lU, number_of_modes);
    lSS = boost::get<0>(res);
    lUU = boost::get<1>(res);
  }

  DoubleMatrix rSS;
  DoubleMatrix rUU;
  if (right_is_enm) {
    DMDuple res = transformENM(rS, rU, number_of_modes);
    rSS = boost::get<0>(res);
    rUU = boost::get<1>(res);
  } else {

    DMDuple res = firstColumns(rS, rU, number_of_modes);
    rSS = boost::get<0>(res);
    rUU = boost::get<1>(res);
  }



  if (square_left)
    for (uint j=0; j<lSS.rows(); ++j)
      lSS[j] *= lSS[j];

  if (square_right)
    for (uint j=0; j<rSS.rows(); ++j)
      rSS[j] *= rSS[j];


  for (uint j=0; j<rSS.rows(); ++j) {
    rSS[j] *= rscale;
    lSS[j] *= lscale;
  }

  if (scale_power)
    rSS = scalePower(lSS, rSS);


  cout << "Covariance Modes: " << number_of_modes << endl;
  if (ntries == 0) {
    double overlap = covarianceOverlap(lSS, lUU, rSS, rUU);
    cout << "Covariance overlap: " << overlap << endl;
  } else {
    if (seed == 0)
      randomSeedRNG();
    else
      rng_singleton().seed(seed);

    boost::tuple<double,double,double> overlap = zCovarianceOverlap(lSS, lUU, rSS, rUU, ntries);
    cout << "Covariance overlap: " << boost::get<1>(overlap) << endl;
    cout << "Z-score: " << boost::get<0>(overlap) << endl;
  }

  double subover = subspaceOverlap(lUU, rUU, subspace_size);
  cout << "Subspace Modes: " << subspace_size << endl;
  cout << "Subspace overlap: " << subover << endl;

}
