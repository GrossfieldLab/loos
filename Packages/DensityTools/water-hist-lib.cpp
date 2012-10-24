// -------------------------------------------------
// Water Histogram Library
// -------------------------------------------------

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009 Tod D. Romo, Alan Grossfield
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester

  This package (LOOS) is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation under version 3 of the License.

  This package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/



#include <water-hist-lib.hpp>


namespace loos {
  namespace DensityTools {


    void ZClipEstimator::reinitialize(pTraj& traj, const std::vector<uint>& frames) {
        std::vector<GCoord> bdd = getBounds(traj, water_, frames);

        bdd[0] -= 1;
        if (bdd[0][2] > zclip_)
          throw(std::runtime_error("Bulk zclip is too small"));
        bdd[0][2] = zclip_;
        bdd[1] += 1;
    
        GCoord gridsize = bdd[1] - bdd[0] + 1;
        gridsize /= gridres_;
        DensityGridpoint dims;
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


    // ------------------------------------------------------------------------

    // Note: no checks on whether the z-slice is sensible...

    void ZSliceEstimator::reinitialize(pTraj& traj, const std::vector<uint>& frames) {
        std::vector<GCoord> bdd = getBounds(traj, water_, frames);

        bdd[0] -= 1;
        bdd[0][2] = zmin_;
        bdd[1] += 1;
        bdd[1][2] = zmax_;
    
        GCoord gridsize = bdd[1] - bdd[0] + 1;
        gridsize /= gridres_;
        DensityGridpoint dims;
        for (int i=0; i<3; ++i)
          dims[i] = static_cast<int>(floor(gridsize[i] + 0.5));
    
        thegrid.resize(bdd[0], bdd[1], dims);
      }

    void ZSliceEstimator::operator()(const double density) {
        for (AtomicGroup::const_iterator i = water_.begin(); i != water_.end(); ++i) {
          GCoord c = (*i)->coords();
          if (c.z() >= zmin_ && c.z() < zmax_)
            thegrid((*i)->coords()) += density;
        }
      }

    double ZSliceEstimator::bulkDensity(void) const {
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

      double ZSliceEstimator::stdDev(const double mean) const {
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


    // ------------------------------------------------------------------------


    void WaterHistogrammer::setGrid(const GCoord& min, const GCoord& max, const double resolution) {
      GCoord gridsize = max - min + 1;
      gridsize /= resolution;
      DensityGridpoint dims;
      for (int i=0; i<3; ++i)
        dims[i] = static_cast<int>(floor(gridsize[i] + 0.5));
      
      grid_.resize(min, max, dims);
    }


    void WaterHistogrammer::setGrid(pTraj& traj, const std::vector<uint>& frames, const double resolution, const double pad) {

      std::vector<GCoord> bdd(2);
      for (uint i=0; i<frames.size(); ++i) {
        traj->readFrame(i);
        traj->updateGroupCoords(protein_);
        std::vector<GCoord> fbdd = the_filter->boundingBox(protein_);
        for (uint j=0; j<3; ++j) {
          if (fbdd[0][j] < bdd[0][j])
            bdd[0][j] = fbdd[0][j];
          if (fbdd[1][j] > bdd[1][j])
            bdd[1][j] = fbdd[1][j];
        }
      }


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
