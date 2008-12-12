#include <iostream>
#include <string>
#include <iterator>
#include <stdexcept>

#include <boost/format.hpp>


#include "Matrix.hpp"

using namespace std;
using namespace lab;

// Quick formatted output of a matrix...
template<class T>
void show(const T& M, const string& s, const string& fmt) {
  cout << s << endl;
  int m = M.rows();
  int n = M.cols();

  for (int j=0; j<m; ++j) {
    for (int i=0; i<n; ++i)
      cout << boost::format(fmt) % M(j, i);
    cout << endl;
  }
}

// Formatted output of the matrix as a linear array...
template<class T>
void showLinear(const T& M, const string& s, const string& fmt) {
  cout << s << endl;
  long n = M.size();
  for (long i=0; i<n; ++i) {
    cout << boost::format(fmt) % M[i];
    cout << endl;
  }
}



int main(int argc, char *argv[]) {

  // How to write/read a matrix...
  if (argc == 2) {
    string name(argv[1]);
    name += ".asc";
    Matrix<float> M(4,4);
    
    int k=0;
    for (int j=0; j<4; ++j)
      for (int i=0; i<4; ++i)
        M(j,i) = k++;
    
    show(M, "M", "%8.2f");
    
    writeAsciiMatrix(name, M, "Testing");
    
    // First form for reading...
    Matrix<float, ColMajor, SharedArray> A = readAsciiMatrix<float, ColMajor, SharedArray>(name);
    cout << boost::format("Read in a %d x %d matrix.\n") % A.rows() % A.cols();
    show(A, "A", "%8.2f");

    // Second form for reading...
    // Relies on default polices...
    Matrix<float> B;
    readAsciiMatrix(name, B);
    cout << boost::format("(2nd form) read in a %d x %d matrix.\n") % B.rows() % B.cols();
    show(B, "B", "%8.2f");
  }

  // Create a row-major matrix with floating point elements...
  Matrix<float, RowMajor> M(4,4);
  int k=0;
  for (int j=0; j<4; ++j)
    for (int i=0; i<4; ++i, ++k)
      M(j,i) = k;

  showLinear(M, "M (row-major)", "%8.2f");

  // You can also output a matrix in Matlab/Octave m-script syntax
  // using the << operator...

  cout << "M = " << M << endl;


  // Create a col-major matrix with integer elements
  Matrix<int, ColMajor> N(4, 4);
  for (int j=0,k=0; j<4; ++j)
    for (int i=0; i<4; ++i, ++k)
      N(j, i) = k;
  show(N, "N (col-major)", "%8.2f");
  showLinear(N, "N (col-major)", "%8.2f");

   // Reinterpret as row-major
  Matrix<float, ColMajor> MM = reinterpretOrder(M);
  show(MM, "N as row-major", "%8.2f");
  showLinear(MM, "MM (linealy)", "%8.2f");

  // Create a triangular matrix...
  Matrix<int, Triangular> T(4,4);
  for (int j=0,k=0; j<4; ++j)
    for (int i=0; i<=j; ++i, ++k)
      T(j, i) = k;
  show(T, "Triangular", "%4d");

  // Accessing the low-level data
  int *ptr = T.get();
  cout << "Low-level access to T\n";
  long s = T.size();
  for (long i=0; i<s; ++i)
    cout << "\t" << ptr[i] << endl;

  // Test of index range-checking...
  bool exception_caught = false;
  try {
    int v = T(4,4);
  }
  catch(out_of_range& e) {
    exception_caught = true;
  }

  if (!exception_caught)
    cout << "***WARNING***\nWe didn't catch an expected out_of_range exception.\n";

  cout << "* iterator test *\n";
  cout << "T = ";
  copy(T.begin(), T.end(), ostream_iterator<int>(cout, ","));
  cout << endl;

  if (argc == 2) {
    cout << "* Writing Triangular Matrix *\n";
    string name(argv[1]);
    name += ".tri";
    writeAsciiMatrix(name, T, "Testing");
    Matrix<int, Triangular> TT;
    readAsciiMatrix(name, TT);
    cout << boost::format("Read in a %d x %d triangular matrix.\n") % TT.rows() % TT.cols();
    show(TT, "T (from file)", "%8.2f");
  }

  cout << "* Copy test *\n";
  M(1,1) = 3.141;
  Matrix<float,RowMajor> MC = M.copy();
  MC(1,1) = 2.718;
  show(M, "Original (1,1)=pi", "%8.2f");
  show(MC, "Copy (1,1)=e", "%8.2f");

  // Sparse test...
  Matrix<float, RowMajor, SparseArray> S(4,4);
  S(1,1) = 1;
  S(2,2) = 3;
  S(3,3) = 5;
  S(1,3) = 7;
  cout << "* Sparse test *\n";
  cout << "actualSize = " << S.actualSize() << endl;
  show(S, "Sparse", "%8.2f");
  cout << "actualSize = " << S.actualSize() << endl;

  // This tests writing/reading of sparse matrices...
  if (argc == 2) {
    string name(argv[1]);
    name += ".spm";

    cout << "* Sparse IO test *\n";
    writeAsciiMatrix(name, S, "Testing");
    Matrix<float, RowMajor, SparseArray> SS;
    readAsciiMatrix(name, SS);

    cout << boost::format("Read in a %d x %d sparse matrix.\n") % SS.rows() % SS.cols();
    cout << "actualSize = " << SS.actualSize() << endl;
    show(SS, "Sparse (read in from file)", "%8.2f");
  }

  // Test of copyMatrix...
  Matrix<float, RowMajor, SparseArray> SC;
  copyMatrix(SC, S);
  cout << "actualSize of copy = " << SC.actualSize() << endl;
  cout << "Density is " << SC.density() << endl;
  
}
