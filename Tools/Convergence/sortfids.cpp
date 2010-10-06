/*
  sortfids
  (c) 2010 Tod D. Romo, Grossfield Lab, URMC
  
  Structural histogram, a la Lyman &
  Zuckerman, Biophys J (2006) 91:164-172

  Usage- sortfids model selection fids hist newfidname

*/



#include <loos.hpp>

using namespace std;
using namespace loos;


typedef vector<AtomicGroup>                   vGroup;
typedef Math::Matrix<double, Math::RowMajor>  Matrix;



struct Adapter {
  Adapter(const Matrix& M) : A(M) { }

  uint size() const { return(A.rows()); }
  const double& operator[](const uint i) const {
    return(A(i,2));
  }

  const Matrix& A;

};


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  
  if (argc != 6) {
    cerr << "Usage- sortfids model sel fids hist newfids\n";
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
