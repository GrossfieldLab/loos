/*
  water-lib.hpp

  (c) 2009 Tod D. Romo, Grossfield Lab, URMC


  Library code common to the water suite...

*/




#if !defined(WATER_LIB_HPP)
#define WATER_LIB_HPP


#include <loos.hpp>


namespace loos {
  namespace DensityTools {

    //! Get the max bounding box for a group over the trajectory
    std::vector<GCoord> getBounds(pTraj& traj, AtomicGroup& group, const std::vector<uint>& frames);

  };

};





#endif
