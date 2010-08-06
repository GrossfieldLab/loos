//
// vsa-lib
//
// (c) 2010 Tod D. Romo, Grossfield Lab, URMC
//


#if !defined(VSA_LIB_HPP)
#define VSA_LIB_HPP

#if defined(__linux__)
extern "C" {
  void dsygvx_(int*, char*, char*, char*, int*, double*, int*, double*, int*, double*, double*, int*, int*, double*, int*, double*, double*, int*, double*, int*, int*, int*, int*);
  void dpotrf_(char*, int*, double*, int*, int*);
  void dtrmm_(char*, char*, char*, char*, int*, int*, double*, double*, int*, double*, int*);
}
#endif




class VSA : public ElasticNetworkModel {
public:
  VSA(SuperBlock* blocker, const uint subn) : 
    ElasticNetworkModel(blocker),
    subset_size_(subn)
  {}

  VSA(SuperBlock* blocker, const uint subn, const loos::DoubleMatrix& M) :
    ElasticNetworkModel(blocker),
    subset_size(subn),
    masses_(M)
  {}


  void solve() {
    buildHessian();
    
    uint n = hessian_.cols();
    uint l = subset_size_ * 3;
    
    DoubleMatrix Hss = submatrix(hessian_, Range(0,l), Range(0,l));
    DoubleMatrix Hee = submatrix(hessian_, Range(l, n), Range(l, n));
    DoubleMatrix Hse = submatrix(hessian_, Range(0,l), Range(l, n));
    DoubleMatrix Hes = submatrix(hessian_, Range(l, n), Range(0, l));

    DoubleMatrix Heei = Math::invert(Hee);
  
    // Build the effective Hessian
    Hssp_ = Hss - Hse * Heei * Hes;

    // Shunt in the event of using unit masses...  We can use the SVD to
    // to get the eigenpairs from Hssp
    if (masses_.rows() == 0) {
      boost::tuple<DoubleMatrix, DoubleMatrix, DoubleMatrix> svdresult = svd(Hssp_);

      eigenvecs_ = boost::get<0>(svdresult);
      eigenvals_ = boost::get<1>(svdresult);

      reverseColumns(eigenvecs_);
      reverseRows(eigenvals_);
      return;

    }


    // Build the effective mass matrix
    DoubleMatrix Ms = submatrix(masses_, Range(0, l), Range(0, l));
    DoubleMatrix Me = submatrix(masses_, Range(l, n), Range(l, n));

    Msp_ = Ms + Hse * Heei * Me * Heei * Hes;

    if (debug) {
      writeAsciiMatrix(prefix + "_Ms.asc", Ms, hdr, false, sp);
      writeAsciiMatrix(prefix + "_Me.asc", Me, hdr, false, sp);
      writeAsciiMatrix(prefix + "_Msp.asc", Msp_, hdr, false, sp);
    }

    // Run the eigen-decomposition...
    if (verbosity > 0) {
      cerr << "Running eigen-decomposition of " << Hssp_.rows() << " x " << Hssp_.cols() << " matrix ...";
      timer.start();
    }
    boost::tuple<DoubleMatrix, DoubleMatrix> eigenpairs;
    eigenpairs = eigenDecomp(Hssp_, Msp_);

    eigenvals_ = boost::get<0>(eigenpairs);
    DoubleMatrix Us = boost::get<1>(eigenpairs);

    // Need to mass-weight the eigenvectors so they're orthogonal in R3...
    if (verbosity > 0)
      cerr << "mass weighting eigenvectors...";

    eigenvecs_ = massWeight(Us, Msp_);

    if (verbosity > 0) {
      timer.stop();
      cerr << "done\n";
      cerr << timer << endl;
    }
  }
  

  void setMasses(const loos::DoubleMatrix& M) {
    masses_ = M;
  };


  // Free up internal storage...
  void free() {
    masses_.reset();
    Msp_.reset();
    Hssp_.reset();
  }

  
private:
  uint subset_size_;
  DoubleMatrix masses_;

  DoubleMatrix Msp_;
  DoubleMatrix Hssp_;
};


#endif
