/*
  Spring Functions from Elastic Network Tools in LOOS
  (c) 2010 Tod D. Romo, Grossfield Lab, URMC
*/



#include "spring_functions.hpp"


// Factory function for SpringFunction's.  In the future, there will
// be additional argument processing here (i.e. constants for the
// spring functions)...
//
// name[,const1[,const2[,...]]]
//
//
// hca,4.0,1,2,3,4
// power,-2.0
//



std::vector<std::string> splitCommaSeparatedList(const std::string& s) {
  // Split of comma-separated string into a vector of strings...

}

SpringFunction* springFactory(const std::string& spring_desc) {

  std::vector<std::string> list = splitCommaSeparatedList(spring_desc);

  if (list[0] == "distance") {
    if (list.size() > 1)
      return(new DistanceWeighted(parseStringAs<double>(list[1])));
    return(new DistanceWeighted);
  }

  if (spring_desc == "weighted")
    return(new DistanceWeighted);
  if (spring_desc == "exponential")
    return(new ExponentialDistance);
  if (spring_desc == "hca" || spring_desc == "HCA")
    return(new HCA);
  if (spring_desc = "membrane")
    


  std::stringstream oss;
  oss << "Error- unknown spring function '" << spring_desc << "'" << std::endl;
  oss << "Try: \"distance\", \"hca\", \"weighted\", or \"exponential\" " << std::endl;
  throw(std::runtime_error(oss.str()));

}
