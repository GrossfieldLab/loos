/*
  Fit aggregator


  (c) 2010 Tod D. Romo, Grossfield Lab, URMC
*/




#if !defined(FIT_AGGREGATOR)
#define FIT_AGGREGATOR


#include "fitter.hpp"

class FitAggregator {
public:
  void push_back(ENMFitter* p) { fitters.push_back(p); }
  double operator()(const std::vector<double>& v) {
    double sum = 0.0;

    for (std::vector<ENMFitter*>::iterator i = fitters.begin(); i != fitters.end(); ++i)
      sum += (**i)(v);

    sum /= fitters.size();

    std::cout << "* Joint = " << -sum << "\n";
    return(sum);
  }


private:
  std::vector<ENMFitter*> fitters;


};



#endif
