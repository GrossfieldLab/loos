#include <GridUtils.hpp>


using namespace std;



namespace loos {
  namespace DensityTools {

    
    std::vector<double> gaussian1d(const int w, const double sigma) {
    
      double a = 1.0 / sqrt(2.0 * M_PI * sigma);
      double b = -1.0 / (2.0 * sigma);
    
      std::vector<double> kernel;
      for (int i=0; i<=w; ++i) {
        double x = 2.0*i/w-1.0;
        kernel.push_back(a*exp(b*x*x));
      }
    
      return(kernel);
    }


  };
};



