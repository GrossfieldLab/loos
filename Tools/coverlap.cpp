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
bool scale_to_svals;
bool scale_to_sum;
bool squares;
uint number_of_modes;
double lscale;
double rscale;
uint subspace_size;

uint skip;


void parseArgs(int argc, char *argv[]) {

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("skip,i", po::value<uint>(&skip)->default_value(6), "# of eigenvalues to skip for ENM")
      ("left_enm,e", po::value<bool>(&left_is_enm)->default_value(false), "Left side contains ENM results")
      ("right_enm,E", po::value<bool>(&right_is_enm)->default_value(false), "Right side contains ENM results")
      ("square_left,s", po::value<bool>(&square_left)->default_value(false), "Square left side (assumes PCA)")
      ("square_right,S", po::value<bool>(&square_right)->default_value(false), "Square right side (assumes PCA)")
      ("scale,r", po::value<bool>(&scale_to_svals)->default_value(false), "Scale ENM eigenvalues (right) to PCA svals (left)")
      ("sum,R", po::value<bool>(&scale_to_sum)->default_value(false), "Scale ENM eigenvalues (right) to PCA svals (left) using sum")
      ("squares,q", po::value<bool>(&squares)->default_value(false), "Use square in sum")
      ("modes,m", po::value<uint>(&number_of_modes)->default_value(0), "Number of modes to compare...  0 = all")
      ("left_scale,k", po::value<double>(&lscale)->default_value(1.0), "Scale left eigenvalues by this constant")
      ("right_scale,K", po::value<double>(&rscale)->default_value(1.0), "Scale right eigenvalues by this constant")
      ("subspace,u", po::value<uint>(&subspace_size)->default_value(0), "# of modes to use for the subspace overlap (0 = same as covariance)");


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

    if (vm.count("help") || !(vm.count("ls") && vm.count("lu") && vm.count("rs") && vm.count("ru"))) {
      cerr << "Usage- " << argv[0] << " [options] ls lU rs rU >output\n";
      cerr << generic;
      exit(-1);
    }

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}

typedef boost::tuple<RealMatrix, RealMatrix>   RMDuple;


RMDuple transformENM(const RealMatrix& S, const RealMatrix& U, const uint nmodes) {
  RealMatrix SS(nmodes, 1);
  RealMatrix UU(U.rows(), nmodes);

  for (uint i=skip; i<nmodes+skip; ++i) {
    SS[i-skip] = 1.0 / S[i];
    for (uint j=0; j<U.rows(); ++j)
      UU(j,i-skip) = U(j,i);
  }

  return(RMDuple(SS,UU));
}


RMDuple firstColumns(const RealMatrix& S, const RealMatrix& U, const uint nmodes) {
  RealMatrix SS(nmodes, 1);
  RealMatrix UU(U.rows(), nmodes);

  for (uint i=0; i<nmodes; ++i) {
    SS[i] = i < S.rows() ? S[i] : 0.0;

    for (uint j=0; j<U.rows(); ++j)
      UU(j, i) = U(j, i);
  }

  return(RMDuple(SS,UU));
}


RealMatrix scaleSvals(const RealMatrix& A, const RealMatrix& B) {
  RealMatrix E(A.rows(), 1);

  double mean = 0.0;
  for (uint j=0; j<A.rows(); ++j)
    mean += B[j] / A[j];

  mean /= A.rows();
  cerr << "Scale factor " << 1/mean << endl;
  for (uint j=0; j<A.rows(); ++j)
    E[j] = B[j] / mean;

  return(E);
}


RealMatrix scaleSquares(const RealMatrix& A, const RealMatrix& B) {

  double sumB = 0.0;
  double sumA = 0.0; 
  for (uint j=0; j<B.rows(); ++j) {
    sumB += B[j];
    if (squares)
      sumA += A[j] * A[j];
    else
      sumA += A[j];
  }

  double scale = sumA / sumB;
  cerr << "Scale factor = " << scale << endl;
  RealMatrix E(B.rows(), 1);
  for (uint j=0; j<B.rows(); ++j)
    E[j] = B[j] * scale;

  return(E);
}




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  parseArgs(argc, argv);

  cerr << "Reading left side matrices...\n";
  RealMatrix lS;
  readAsciiMatrix(lefts_name, lS);
  RealMatrix lU;
  readAsciiMatrix(leftU_name, lU);
  cerr << boost::format("Read in %d x %d eigenvectors...\n") % lU.rows() % lU.cols();
  cerr << boost::format("Read in %d eigenvalues...\n") % lS.rows();

  cerr << "Reading in right side matrices...\n";
  RealMatrix rS;
  readAsciiMatrix(rights_name, rS);
  RealMatrix rU;
  readAsciiMatrix(rightU_name, rU);
  cerr << boost::format("Read in %d x %d eigenvectors...\n") % rU.rows() % rU.cols();
  cerr << boost::format("Read in %d eigenvalues...\n") % rS.rows();

  if (number_of_modes == 0) {
    number_of_modes = lS.rows() > rS.rows() ? lS.rows() : rS.rows();
    if (left_is_enm || right_is_enm)
      number_of_modes -= skip;
  }


  if (subspace_size == 0)
    subspace_size = number_of_modes;
  else
    if (subspace_size > number_of_modes) {
      cerr << "ERROR- subspace size cannot exceed number of modes for covariance overlap\n";
      exit(-1);
    }

  RealMatrix lSS;
  RealMatrix lUU;
  if (left_is_enm) {
    RMDuple res = transformENM(lS, lU, number_of_modes);
    lSS = boost::get<0>(res);
    lUU = boost::get<1>(res);
  } else {
    RMDuple res = firstColumns(lS, lU, number_of_modes);
    lSS = boost::get<0>(res);
    lUU = boost::get<1>(res);
  }

  RealMatrix rSS;
  RealMatrix rUU;
  if (right_is_enm) {
    RMDuple res = transformENM(rS, rU, number_of_modes);
    rSS = boost::get<0>(res);
    rUU = boost::get<1>(res);
  } else {

    RMDuple res = firstColumns(rS, rU, number_of_modes);
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

  if (scale_to_svals)
    rSS = scaleSvals(lSS, rSS);
  else if (scale_to_sum)
    rSS = scaleSquares(lSS, rSS);

  double overlap = covarianceOverlap(lSS, lUU, rSS, rUU);
  double subover = subspaceOverlap(lUU, rUU, subspace_size);

  cout << "Modes: " << number_of_modes << endl;
  cout << "Covariance overlap: " << overlap << endl;
  cout << "Modes: " << subspace_size << endl;
  cout << "Subspace overlap: " << subover << endl;

}
