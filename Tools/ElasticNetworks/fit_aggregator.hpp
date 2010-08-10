/*
  Fit aggregator


  (c) 2010 Tod D. Romo, Grossfield Lab, URMC
*/




#if !defined(FIT_AGGREGATOR)
#define FIT_AGGREGATOR


#include "fitter.hpp"

class FitAggregator {
public:
  FitAggregator() : iters_(0), verbose_(true) { }

  bool verbose() const { return(verbose_); }
  void verbose(const bool b) { verbose_ = b; }

  uint iterations() const { return(iters_); }


  void push_back(ENMFitter* p) { fitters.push_back(p); }

  double operator()(const std::vector<double>& v) {
    double sum = 0.0;

    for (std::vector<ENMFitter*>::iterator i = fitters.begin(); i != fitters.end(); ++i)
      sum += (**i)(v);

    sum /= fitters.size();
    
    ++iters_;
    if (verbose_) 
      std::cout << "* (" << iters_ << ") Joint = " << -sum << "\n";

    return(sum);
  }

  void resetCount() { iters_ = 0; }


private:
  uint iters_;
  bool verbose_;
  std::vector<ENMFitter*> fitters;

};



#endif
