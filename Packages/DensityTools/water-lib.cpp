/*
  water-lib.cpp

  (c) 2009 Tod D. Romo, Grossfield Lab, URMC


  Library code common to the water suite...

*/




#include "water-lib.hpp"


namespace banal {
  namespace water {

    std::vector<GCoord> getBounds(pTraj& traj, AtomicGroup& g, const std::vector<uint>& indices) {
      double d = std::numeric_limits<double>::max();
      GCoord min(d, d, d);
      GCoord max(-d, -d, -d);
      
      for (std::vector<uint>::const_iterator i = indices.begin(); i != indices.end(); ++i) {
        traj->readFrame(*i);
        traj->updateGroupCoords(g);
        std::vector<GCoord> bdd = g.boundingBox();
        for (int j=0; j<3; ++j) {
          if (bdd[0][j] < min[j])
            min[j] = bdd[0][j];
          if (bdd[1][j] > max[j])
            max[j] = bdd[1][j];
        }
      }
      
      std::vector<GCoord> bdd;
      bdd.push_back(min);
      bdd.push_back(max);
      return(bdd);
    }
 

  };

};
