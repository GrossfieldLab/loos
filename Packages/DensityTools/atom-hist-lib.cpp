/*
  water-hist-lib.cpp

  
  (c) 2009 Tod D. Romo, Grossfield Lab, URMC

  Water histogram library
*/




#include "water-hist-lib.hpp"


namespace banal {
  namespace water {


    void ZClipEstimator::reinitialize(pTraj& traj, const std::vector<uint>& frames) {
        std::vector<GCoord> bdd = banal::water::getBounds(traj, water_, frames);

        bdd[0] -= 1;
        if (bdd[0][2] > zclip_)
          throw(std::runtime_error("Bulk zclip is too small"));
        bdd[0][2] = zclip_;
        bdd[1] += 1;
    
        GCoord gridsize = bdd[1] - bdd[0] + 1;
        gridsize /= gridres_;
        lab::SGridpoint dims;
        for (int i=0; i<3; ++i)
          dims[i] = static_cast<int>(floor(gridsize[i] + 0.5));
    
        thegrid.resize(bdd[0], bdd[1], dims);
      }

    void ZClipEstimator::operator()(const double density) {
        for (AtomicGroup::const_iterator i = water_.begin(); i != water_.end(); ++i) {
          GCoord c = (*i)->coords();
          if (c.z() >= zclip_)
            thegrid((*i)->coords()) += density;
        }
      }

    double ZClipEstimator::bulkDensity(void) const {
        double mean = 0.0;
        long n = 0;
        for (long i = 0; i < thegrid.maxGridIndex(); ++i) {
          if (!count_zero && thegrid(i) == 0.0)
            continue;
          mean += thegrid(i);
          ++n;
        }

        return(mean/n);
      }

      double ZClipEstimator::stdDev(const double mean) const {
        double std = 0.0;
        long n = 0;
        for (long i = 0; i < thegrid.maxGridIndex(); ++i)
          if (thegrid(i) != 0.0) {
            double d = thegrid(i) - mean;
            std += d*d;
            ++n;
          }

        std = sqrt(std/(n-1));
        return(std);
      }


    void WaterHistogrammer::setGrid(GCoord min, GCoord max, const double resolution) {
      GCoord gridsize = max - min + 1;
      gridsize /= resolution;
      lab::SGridpoint dims;
      for (int i=0; i<3; ++i)
        dims[i] = static_cast<int>(floor(gridsize[i] + 0.5));
      
      grid_.resize(min, max, dims);
    }


    void WaterHistogrammer::setGrid(pTraj& traj, const std::vector<uint>& frames, const double resolution, const double pad) {
      std::vector<GCoord> bdd = banal::water::getBounds(traj, protein_, frames);
      setGrid(bdd[0] - pad, bdd[1] + pad, resolution);
    }
    
    
      void WaterHistogrammer::accumulate(const double density) {
        std::vector<int> picks = the_filter->filter(water_, protein_);
        for (uint i = 0; i<picks.size(); ++i)
          if (picks[i]) {
            GCoord c = water_[i]->coords();

            if (!grid_.inRange(grid_.gridpoint(c)))
              ++out_of_bounds;
            else
              grid_(c) += density;
          }
        (*estimator_)(density);
      }

      void WaterHistogrammer::accumulate(pTraj& traj, const std::vector<uint>& frames) {
        estimator_->reinitialize(traj, frames);
        double density = 1.0 / frames.size();
        for (std::vector<uint>::const_iterator i = frames.begin(); i != frames.end(); ++i) {
          traj->readFrame(*i);
          traj->updateGroupCoords(protein_);
          traj->updateGroupCoords(water_);

          accumulate(density);
        }
      }

    };
};
