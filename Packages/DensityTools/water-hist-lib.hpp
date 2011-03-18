/*
  water-hist-lib.hpp

  (c) 2009 Tod D. Romo, Grossfield Lab, URMC

  Water histogram library
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




    class WaterHistogrammer {
    public:
      WaterHistogrammer(const AtomicGroup& protein, const AtomicGroup& water, BulkEstimator* est, WaterFilter::Base* filter) :
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
      WaterFilter::Base* the_filter;
      long out_of_bounds;
      DensityGrid<double> grid_;
    };


  };

};



#endif
