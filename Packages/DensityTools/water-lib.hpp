/*
  water-lib.hpp

  (c) 2009 Tod D. Romo, Grossfield Lab, URMC


  Library code common to the water suite...

*/




#if !defined(WATER_LIB_HPP)
#define WATER_LIB_HPP


#include <loos.hpp>


namespace banal {
  namespace water {

    using namespace loos;

    std::vector<GCoord> getBounds(pTraj&, AtomicGroup&, const std::vector<uint>&);

  };

};





#endif
