#if !defined(MC_FIT_HPP)
#define MC_FIT_HPP

#include <iostream>
#include <vector>
#include <ctime>            
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

typedef boost::minstd_rand rand_num;

template<typename T = double>
class mcoptimo {
  typedef std::vector< std::vector<T> >     VVectors;

public:
   void reseed(const &previous, &current){
    current.setParams() = randomize(previous.getParams());
  }

  void accept(&previous, const &current){
    if (current.getCov() > previous.getCov()){       //accept
      previous = current;
      // }else if (current.getCov() > e^(Energy/T_star)){ //accept
      // previous = current;
    } // else{                                           //reject
      //   continue; 
      // }
  }

  void optimum(const &current, &best){
    if (current.getCov() > best.getCov()) best = current;
  }

  void randomize(const &previous, &current){
    //use the BOOST random generator per alan's suggestion
    //using type pseudo-random number generator, based on the previous k's
    rand_num generator(time(0)); //<---Talk to Tod about other seeds!!!
    vector<double> newParams;
    for (int i = 0; i < previous.getParams->paramSize(); ++i){
      double springval = 0;//add stuff here to grab that a particular param 
      double lowerbound = springval / 2;
      double upperbound = 3 * springval / 2;
      boost::uniform_real<> uni_dist(lowerbound , upperbound);
      boost::variate_generator<rand_num&, boost::uniform_real<> > uni(generator, uni_dist);
      newParams.push_back(uni() + uni());
    }
    current.setParams(newParams);
  }

  double getCov(){return this->cov;}
  void setCov(&cov_){this->cov = cov_}

  params getParams(){return this->K;}
  void setParams(&K_/*a vector of params!*/){
    for (int i = 0; i < K_.size(); ++i)
      K.push_back(K_);
  }

  void operator=(&A, &B){
    A.setCov(B.getCov);
    A.setParams(B.getParams);
  }
  void operator()(){}//i have no idea why i need this...but all ftors use it!!!

private:
  double cov;
  params K; //should this be a springFactory instance???



  //////////////////////////////////////////////////



  // The core of the optimizer...
  template<class C>
  void core(C& ftor) {
    int i, next_worst, j, mpts;
    T best, current, previous;//sum, saved, val, num, den;

    mpts = ndim + 1;

    for (j=0; j<ndim; j++) {
      for (sum = 0.0, i = 0; i<mpts; i++)

        sum += simplex[i][j];
      simpsum[j] = sum;
    }

    int n_evals = 0;

    while (n_evals <= maxiters) {
      best = 1;
      if (values[0] > values[1]) {
        worst = 0;
        next_worst = 1;
      } else {
        worst = 1;
        next_worst = 0;
      }


      // Find candidate vertices in the simplex...

      for (i=0; i<mpts; i++) {
        if (values[i] <= values[best])
          best = i;
        if (values[i] > values[worst]) {
          next_worst = worst;
          worst = i;
        } else if (values[i] > values[next_worst] && i != worst)
          next_worst = i;
      }

      // Check for convergence...
      // Some conditions may have the numerator and denominator equal (or both 0)
      // which can cause problems below, hence the special test.

      num = fabs(values[worst] - values[best]);
      den = fabs(values[worst]) + fabs(values[best]);
      rtol = 2.0 * num / den;
      if (rtol < tol || num == den)
        return;

      // Now try reflecting, contracting, expanding the simplex...

      n_evals += 2;
      val = modify(-1.0, ftor);
      if (val <= values[best])
        val = modify(2.0, ftor);
      else if (val >= values[next_worst]) {

        saved = values[worst];
        val = modify(0.5, ftor);

        if (val >= saved) {
          for (i=0; i<mpts; i++) {
            if (i != best) {
              for (j=0; j<ndim; j++)
                simplex[i][j] = simpsum[j] = 0.5 * (simplex[i][j] + simplex[best][j]);
              values[i] = ftor(simpsum);
            }
          }
          n_evals += ndim;

          for (j=0; j<ndim; j++) {
            for (sum=0.0, i=0; i<mpts; i++)
              sum += simplex[i][j];
            simpsum[j] = sum;
          }
        }

      } else
        --n_evals;

    }

  }


  void allocateSpace(const int n) {
    q = std::vector<T>(n+1, 0.0);
    qq = std::vector<T>(n+1, 0.0);
    simpsum = std::vector<T>(n, 0.0);
    values = std::vector<T>(n+1, 0.0);
    trial = std::vector<T>(n, 0.0);
    simplex.clear();
    for (int i=0; i<n+1; ++i)
      simplex.push_back(std::vector<T>(n, 0.0));
  }



public:

  Simplex(const int n) : tol(1e-3), ndim(n), maxiters(2000), best(-1), worst(-1) {
    allocateSpace(n);
  }

  //! Set the number of dimensions
  void dim(const int n) { ndim = n; allocateSpace(n); }

  //! Initial guess
  void seedLengths(const std::vector<T>& seeds) { characteristics = seeds; }

  //! Convergence criterion
  void tolerance(const double d) { tol = d; }

  //! Limit on the number of function evaluations to perform
  void maximumIterations(const int n) { maxiters = n; }
  
  //! Retrieve the final (best fit) parameters
  std::vector<T> finalParameters(void) const {
    if (best < 0)
      throw(std::logic_error("Simplex has not been optimized"));
    return(simplex[best]);
  }

  //! Final (best) value
  T finalValue(void) const { return(values[best]); }
 
  //! Optimize via a functor
  template<class C>
  std::vector<T> optimize(std::vector<T>& f, C& ftor) {
    int i, j, n;

    if (characteristics.size() != f.size())
      throw(std::logic_error("Invalid seed"));

    n = ndim + 1;

    // Initial simplex is based on Nelder-Meade's construction...

    for (i=0; i<ndim; i++) {
      q[i] = characteristics[i] * ( (sqrtf(n) + ndim) / (n * sqrt(2.0)) );
      qq[i] = characteristics[i] * ( (sqrtf(n) - 1.0) / (n * sqrt(2.0)) );
    }

    for (j=0; j<n; j++)
      for (i=0; i<ndim; i++)
        if (j == i+1)
          simplex[j][i] = f[i] + qq[i];
        else
          simplex[j][i] = f[i] + q[i];

    for (j=0; j<n; ++j) {
      values[j] = ftor(simplex[j]);
    }

    core(ftor);
    return(simplex[best]);
  }

private:
  double tol;
  int ndim, maxiters, best, worst;
  double rtol;

  std::vector<T> characteristics;
  std::vector<T> simpsum;
  std::vector<T> values;
  std::vector<T> q;
  std::vector<T> qq;
  std::vector<T> trial;
  VVectors simplex;

};




#endif
