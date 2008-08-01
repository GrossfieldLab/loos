#include <loos.hpp>
#include <boost/format.hpp>




GMatrix randomMatrix(const double range) {
  GMatrix M;

  loos::base_generator_type& rng = loos::rng_singleton();
  boost::uniform_real<> rngmap(-range, range);
  boost::variate_generator<loos::base_generator_type&, boost::uniform_real<> > randoms(rng, rngmap);

  for (int i=0; i<16; i++) {
    M[i] = randoms();
  }

  return(M);
}


GMatrix manualMultiply(const GMatrix& A, const GMatrix& B) {
  GMatrix C;
  
  for (int j=0; j<4; j++) {
    for (int i=0; i<4; i++) {
      double d = 0;
      for (int k=0; k<4; k++) {
	d += A(j,k) * B(k,i);
      }
      C(j,i) = d;
    }
  }

  return(C);
}

double rmsd(const GMatrix& A, const GMatrix& B) {
  double m = 0.0, msd = 0.0;
  GMatrix C;
  
  for (int i=0; i<16; i++)
    C[i] = A[i] - B[i];

  for (int i=0; i<16; i++)
    m += C[i];
  m /= 16;

  for (int i=0; i<16; i++)
    msd += (C[i] - m) * (C[i] - m);
  
  msd /= 16;
  return(sqrt(msd));
  
}


int main(int argc, char *argv[]) {

  long niters = atol(argv[1]);
  double range = strtod(argv[2], 0);
  double tol = strtod(argv[3], 0);
  double avg = 0.0;

  GMatrix C;
  GMatrix C2;

  for (long i=0; i<niters; i++) {
    GMatrix A = randomMatrix(range);
    C *= A;
    GMatrix M;
    M = manualMultiply(C2, A);
    C2 = M;

    double err = rmsd(C, C2);
    //    cerr << "Probe: " << err << endl;
    avg += err;
    if (isnan(err) || err >= tol) {
      cerr << "Failuer at iteration " << i << endl;
      if (isnan(err))
	cerr << "NaN returned from rmsd!\n";
      else
	cerr << boost::format("Failure with RMSD=%lf for the following matrices:\n") % err;
      cerr << "A:\n" << A << endl;
      cerr << "C:\n" << C << endl;
      cerr << "C2:\n" << C2 << endl;
      exit(-99);
    }
  }

  cout << boost::format("%d iterations with tol %lf passed.\n") % niters % tol;
  cout << boost::format("Sum of rmsd = %lf\n") % avg;
  cout << boost::format("Average rmsd was %lf\n") % (avg / static_cast<double>(niters));

}


