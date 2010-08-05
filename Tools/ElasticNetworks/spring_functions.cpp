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
// hca,4.0,1,2,3,4
// power,-2.0
//
void splitCommaSeparatedList(const std::string& s, std::vector<std::string>& holder){
  std::string::size_type prev =s.find_first_not_of(",", 0);
  std::string::size_type pos = s.find_first_of(",", prev);
    
  while (pos != std::string::npos || prev != std::string::npos){
    holder.push_back(s.substr(prev, pos-prev));
    prev = s.find_first_not_of(",", pos);
    pos = s.find_first_of(",", prev);
  } 
}


SpringFunction* springFactory(const std::string& spring_desc) {

  std::vector<std::string> list = splitCommaSeparatedList(spring_desc);

  if (list[0] == "distance") {
    if (list.size() > 1)
      return(new DistanceCutoff(parseStringAs<double>(list[1])));//will this jsut send list[1] or does parseStringAs concat??  this won't work it list[1] is not a double
    return(new DistanceCutoff);
  }

  if (list[0] == "weighted"){
    if (list.size() > 1)
      return(new DistanceWeighted(parseStringAs<double>(list[1])));
    return(new DistanceWeighted);
  }

  if (list[0] == "exponential"){
    //i don't think ExponentialDistance takes any additional args at this time....
    // if (list.size() > 1)
    //   return(new ExponentialDistance(parseStringAs<double>(list[1])));
    return(new ExponentialDistance);
  }
  
  if (list[0] == "hca" || spring_desc == "HCA"){
    if (list.size() > 1)
      return(new HCA(parseStringAs<double>(list[1]),parseStringAs<double>(list[2]),parseStringAs<double>(list[3]),parseStringAs<double>(list[4]),parseStringAs<double>(list[5])));
    return(new HCA);
  }
    


  std::stringstream oss;
  oss << "Error- unknown spring function '" << spring_desc << "'" << std::endl;
  oss << "Try: \"distance\", \"hca\", \"weighted\", or \"exponential\" " << std::endl;
  throw(std::runtime_error(oss.str()));

}
