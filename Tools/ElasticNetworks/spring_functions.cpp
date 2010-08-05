/*
  Spring Functions from Elastic Network Tools in LOOS
  (c) 2010 Tod D. Romo, Grossfield Lab, URMC
*/



#include "spring_functions.hpp"


// Factory function for SpringFunction's.  In the future, there will
// be additional argument processing here (i.e. constants for the
// spring functions)...

SpringFunction* springFactory(const std::string& spring_desc) {

  if (spring_desc == "distance")
    return(new DistanceCutoff);
  if (spring_desc == "weighted")
    return(new DistanceWeighted);
  if (spring_desc == "exponential")
    return(new ExponentialDistance);
  if (spring_desc == "hca" || spring_desc == "HCA")
    return(new HCA);


  std::stringstream oss;
  oss << "Error- unknown spring function '" << spring_desc << "'" << std::endl;
  oss << "Try: \"distance\", \"hca\", \"weighted\", or \"exponential\" " << std::endl;
  throw(std::runtime_error(oss.str()));
}

SpringFunction* springFactory(const std::string& spring_desc, const double a, const double b, const double c, const double d){
  if (spring_desc == "hca" || spring_desc == "HCA" ){


    return(new HCA);
  }

}
