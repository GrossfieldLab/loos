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
std::vector<std::string> splitCommaSeparatedList(const std::string& s){
  std::vector<std::string> holder;
  std::string::size_type prev =s.find_first_not_of(",", 0);
  std::string::size_type pos = s.find_first_of(",", prev);
    
  while (pos != std::string::npos || prev != std::string::npos){
    holder.push_back(s.substr(prev, pos-prev));
    prev = s.find_first_not_of(",", pos);
    pos = s.find_first_of(",", prev);
  } 

  return(holder);
}


SpringFunction* springFactory(const std::string& spring_desc) {

  std::vector<std::string> list = splitCommaSeparatedList(spring_desc);
  SpringFunction *spring;

  if (list[0] == "distance")
    spring = new DistanceCutoff;
  else if (list[0] == "weighted")
    spring = new DistanceWeight;
  else if (list[0] == "exponential")
    spring = new ExponentialDistance;
  else if  (list[0] == "hca" || spring_desc == "HCA")
    spring = new HCA;
  else {
    std::stringstream oss;
    oss << "Error- unknown spring function '" << spring_desc << "'" << std::endl;
    oss << "Try: \"distance\", \"hca\", \"weighted\", or \"exponential\" " << std::endl;
    throw(std::runtime_error(oss.str()));
  }

  if (list.size() > 1) {
    std::vector<double> params;

    for (uint i=1; i <= spring->paramSize(); ++i)
      params.push_back(loos::parseStringAs<double>(list.at(i)));
    spring->setParams(params);
  }

  return(spring);
}
