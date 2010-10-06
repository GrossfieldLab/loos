/*
  Fiducial library

  (c) 2010 Tod D. Romo, Grossfield Lab, URMC
*/


#if !defined(FIDLIB_HPP)
#define FIDLIB_HPP



#include <loos.hpp>


typedef std::vector<int>                   vecInt;
typedef std::vector<uint>                 vecUint;
typedef std::vector<loos::AtomicGroup>   vecGroup;
typedef std::vector<double>             vecDouble;





vecUint findFreeFrames(const vecInt& map);

boost::tuple<vecInt, vecUint, vecGroup, vecDouble> assignFrames(loos::AtomicGroup& model, loos::pTraj& traj, const vecUint& frames, const double f);

int findMaxBin(const vecInt& assignments);
vecUint histogramBins(const vecInt& assignments);



#endif
