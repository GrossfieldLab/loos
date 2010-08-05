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

  std::stringstream oss;
  oss << "Error- unknown spring function '" << spring_desc << "'" << std::endl;
  throw(std::runtime_error(oss.str()));
}

