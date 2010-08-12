/*
  Spring Functions from Elastic Network Tools in LOOS
  (c) 2010 Tod D. Romo, Grossfield Lab, URMC
*/


#include "spring_functions.hpp"


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



/** Factory function for spring constants.  The spring description is
 * the name of the spring function with an optional comma-separated
 * list of parameters to pass to it, i.e.
 *    distance
 *    distance,15.0
 *    hca,1,2,3,4,5
 *
 * Will throw in the event of an error.
 */
SpringFunction* springFactory(const std::string& spring_desc) {

  std::vector<std::string> list = splitCommaSeparatedList(spring_desc);
  SpringFunction *spring = 0;

  if (list[0] == "distance")
    spring = new DistanceCutoff;
  else if (list[0] == "weighted")
    spring = new DistanceWeight;
  else if (list[0] == "exponential")
    spring = new ExponentialDistance;
  else if  (list[0] == "hca" || spring_desc == "HCA")
    spring = new HCA;
  else
    throw(BadSpringFunction("Bad Spring Function Name"));

  if (list.size() > 1) {
    if (list.size() < spring->paramSize())
      throw(BadSpringParameter("Too few spring parameters"));

    std::vector<double> params;

    for (uint i=1; i <= spring->paramSize(); ++i)
      params.push_back(loos::parseStringAs<double>(list.at(i)));

    spring->setParams(params);
    if (!spring->validParams())
      throw(BadSpringParameter("Bad Spring Parameter"));
  }

  return(spring);
}



// Returns a list of names for spring functions.  Needs to match
// what's checked for above...

std::vector<std::string> springNames() {
  std::vector<std::string> names;

  names.push_back("distance");
  names.push_back("weighted");
  names.push_back("exponential");
  names.push_back("hca");

  return(names);
}
