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




#if !defined(WATER_HIST_LIB_HPP)
#define WATER_HIST_LIB_HPP

#include <iostream>

#include <DensityGrid.hpp>
#include <GridUtils.hpp>
#include <water-lib.hpp>
#include <internal-water-filter.hpp>


namespace loos {

  namespace DensityTools {


    class BulkEstimator {
    public:
      virtual ~BulkEstimator() { }
      virtual void reinitialize(pTraj&, const std::vector<uint>&) =0;
      virtual void operator()(const double) =0;
      virtual double bulkDensity(void) const =0;
      virtual double stdDev(const double) const =0;
      virtual void clear() =0;

      friend std::ostream& operator<<(std::ostream& os, const BulkEstimator& b) {
        return(b.print(os));
      }

    private:
      virtual std::ostream& print(std::ostream& os) const =0;
    };


    class NullEstimator : public BulkEstimator {
    public:
      void reinitialize(pTraj& p, const std::vector<uint>& f) { }
      void operator()(const double d) { }
      double bulkDensity(void) const { return(1.0); }
      double stdDev(const double d) const { return(0.0); }
      void clear(void) { }

    private:
      std::ostream& print(std::ostream& os) const {
        os << "No bulk estimate";
        return(os);
      }
    };


    class ZClipEstimator : public BulkEstimator {
    public:
      ZClipEstimator(AtomicGroup& water, pTraj& traj, const std::vector<uint>& frames, const double zclip, const double gridres)
        : water_(water), zclip_(zclip), gridres_(gridres), count_zero(false)
      {
        reinitialize(traj, frames);
      }

      void countZero(const bool flag = true) { count_zero = flag; }

      void reinitialize(pTraj& traj, const std::vector<uint>& frames);
      void operator()(const double density);
      double bulkDensity(void) const;
      double stdDev(const double mean) const;
      void clear(void) { thegrid.clear(); }

    private:
      std::ostream& print(std::ostream& os) const {
        os << boost::format("ZClipEstimator = %s x %s @ %s") % thegrid.minCoord() % thegrid.maxCoord() % thegrid.gridDims();
        return(os);
      }

    private:
      AtomicGroup water_;
      double zclip_, gridres_;
      bool count_zero;
      DensityGrid<double> thegrid;
    };


    class ZSliceEstimator : public BulkEstimator {
    public:
      ZSliceEstimator(AtomicGroup& water, pTraj& traj, const std::vector<uint>& frames, const double zmin, const double zmax, const double gridres)
        : water_(water), zmin_(zmin), zmax_(zmax), gridres_(gridres), count_zero(false)
      {
        reinitialize(traj, frames);
      }

      void countZero(const bool flag = true) { count_zero = flag; }

      void reinitialize(pTraj& traj, const std::vector<uint>& frames);
      void operator()(const double density);
      double bulkDensity(void) const;
      double stdDev(const double mean) const;
      void clear(void) { thegrid.clear(); }

    private:
      std::ostream& print(std::ostream& os) const {
        os << boost::format("ZSliceEstimator = %s x %s @ %s") % thegrid.minCoord() % thegrid.maxCoord() % thegrid.gridDims();
        return(os);
      }

    private:
      AtomicGroup water_;
      double zmin_, zmax_, gridres_;
      bool count_zero;
      DensityGrid<double> thegrid;
    };




    class WaterHistogrammer {
    public:
      WaterHistogrammer(const AtomicGroup& protein, const AtomicGroup& water, BulkEstimator* est, WaterFilterBase* filter) :
        protein_(protein), water_(water), estimator_(est), the_filter(filter), out_of_bounds(0) { }

      void clear() { grid_.clear(); out_of_bounds = 0; }

      void setGrid(GCoord min, GCoord max, const double resolution);
      void setGrid(pTraj& traj, const std::vector<uint>& frames, const double resolution, const double pad = 0.0);

      void accumulate(const double density);
      void accumulate(pTraj& traj, const std::vector<uint>& frames);
      DensityGrid<double> grid() const { return(grid_); }
      long outOfBounds() const { return(out_of_bounds); }



    private:
      AtomicGroup protein_, water_;
      BulkEstimator* estimator_;
      WaterFilterBase* the_filter;
      long out_of_bounds;
      DensityGrid<double> grid_;
    };


  };

};



#endif
