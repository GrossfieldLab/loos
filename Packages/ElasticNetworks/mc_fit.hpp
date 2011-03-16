#if !defined(MC_FIT_HPP)
#define MC_FIT_HPP

#include <iostream>
#include <vector>
#include <ctime>            
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

typedef boost::minstd_rand rand_num;

struct ConstantAcceptor {
  ConstantAcceptor() : val(0.25) { }
  ConstantAcceptor(double d) : val(d) { }
  double operator()(const uint iter) { return(val); }

  double val;
};


struct ExponentialAcceptor {
  ExponentialAcceptor() : k(1.0) { }
  ExponentialAcceptor(const double scale) : k(scale) { }
  double operator()(const uint iter) { return(exp(-k * iter)); }

  double k;
};


template<typename T = double>
class mcoptimo {
  typedef std::vector< std::vector<T> >     VVectors;

public:
  //  void reseed(const mcoptimo& previous, mcoptimo &current){
  //    current.setParams(randomize(previous.getParams()));
  //    //    current.setParams() = randomize(previous.getParams());
  // }
  

  void setParams(const vector<T>& v) {
    myparams = v;
  }
  
  vector<T> getParams() const { return(myparams); }


  //  vector<T>& setParams() { return(myparams); }


  //private:
  //  vector<T> myparams;


//   void accept(&previous, const &current){
//     if (current.getCov() > previous.getCov()){       //accept
//       previous = current;
//       // }else if (current.getCov() > e^(Energy/T_star)){ //accept
//       // previous = current;
//     } // else{                                           //reject
//       //   continue; 
//       // }
//   }

  //  void optimum(const &current, &best){
  //    if (current.getCov() > best.getCov()) best = current;
  //  }

  // Return a perturbed set of parameters
  vector<T> randomize(const vector<T>& current, const vector<T>& sizes) {
    //use the BOOST random generator per alan's suggestion
    //using type pseudo-random number generator, based on the previous k's
    //rand_num generator(time(0)); //<---Talk to Tod about other seeds!!!
    //    base_generator_type& generator = rng_singleton();

    vector<T> newParams;
    for (int i = 0; i < previous.getParams->paramSize(); ++i){//i want the size of the vector myparams
                                                              //for the 'previous' instance of mcoptimize
      T springval = 0;//add stuff here to grab that a particular param 

      // Set range/bounds to the size here...
      // T lowerbound = springval / 2;
      // T upperbound = 3 * springval / 2;
      // boost::uniform_real<> uni_dist(lowerbound , upperbound);
      boost::uniform_real<> uni_dist(-1, 1);
      boost::variate_generator<rand_num&, boost::uniform_real<> > uni(generator, uni_dist);
      double scaled_random = (uni() + uni()) * current[i]; 
      newParams.push_back(scaled_random + current[i]);
    }
    return(newParams);
  }

  //  template<class C>
  T randomize(uint iter){
    base_generator_type& single_random = rng_singleton();
    T upperbound = acceptance(iter);
    boost::uniform_real<> uni_dist(0, upperbound*2);
    boost::variate_generator<rand_num&, boost::uniform_real<> > uni(single_random, uni_dist);
    double my_random = uni();
    return(my_random);
  }

//   template<class C>
//   T acceptance(uint iter, const uint stepSize, const uint absTemp){
//     //maybe use scope to grab stepSize and absTemp
//     //but what scope do they belong in??
                                                                  
//     //this should not be hard wired.  
//     //we want another ftor that pts
//     //to an acceptance function
//     T cutoff = absTemp - (stepSize * iter);
//     return(cutoff);
//  }


  template<class C, class Acceptor = ConstantAcceptor>
  vector<T> takeStep(const vector<T>& current, const vector<T>& sizes, C& ftor, Acceptor& acc, uint iter) {
    
    vector<T> newStep = randomize(current, sizes);
    T prev = ftor(current);
    T val = ftor(newStep);

    if (val < prev){
      return(newStep);
    }elseif (randomize(iter) < acc(iter)){//define both of these!!!
      return(newStep);
    }
    return(current);
  }

  

  // vector<double> takeStep(const vector<double>& current, const vector<double>& sizes, FitAggregator& ftor)
  template<class C, class A = ConstantAcceptor>
  vector<T> optimize(const vector<T>& current, C& ftor, A& ator) {
    
    vector<T> params(current);//<<----do i need this??
    vector<T> sizes(initial_sizes);
    uint iter = 0;
    
    while (!converged) {
      vector<T> params = takeStep(params, sizes, ftor, ator, iter);
      ++iter;

    return(params);
  }

  void setSizes(const vector<T>& s) {
    initial_sizes = s;
  }


private:
  //double cov;
  //params K; //should this be a springFactory instance???
  vector<T> initial_sizes;
  vector<T> myparams;


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
