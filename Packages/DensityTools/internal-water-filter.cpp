/*
  internal_water_filter.cpp

  (c) 2008 Tod D. Romo,
  Grossfield Lab,
  University of Rochester Medical and Dental School
*/



#include <internal-water-filter.hpp>


using namespace std;

namespace loos {
  namespace DensityTools {


    string WaterFilterBox::name(void) const {
      stringstream s;
      s << boost::format("WaterFilterBox(pad=%f)") % pad_;
      return(s.str());
    }
    

    vector<int> WaterFilterBox::filter(const AtomicGroup& solv, const AtomicGroup& prot) {
      bdd_ = boundingBox(prot);
      vector<int> result(solv.size());
  
      AtomicGroup::const_iterator ai;
  
      uint j = 0;
      for (ai = solv.begin(); ai != solv.end(); ++ai) {
        bool flag = true;
        GCoord c = (*ai)->coords();
        for (int i=0; i<3; ++i)
          if (c[i] < bdd_[0][i] || c[i] > bdd_[1][i]) {
            flag = false;
            break;
          }

        result[j++] = flag;
      }
      return(result);
    }


    double WaterFilterBox::volume(void) {
      GCoord v = bdd_[1] - bdd_[0];
      return(v[0] * v[1] * v[2]);
    }

  
    vector<GCoord> WaterFilterBox::boundingBox(const AtomicGroup& grp) {
      vector<GCoord> bdd = grp.boundingBox();
      bdd[0] = bdd[0] - pad_;
      bdd[1] = bdd[1] + pad_;

      return(bdd);
    }


    // --------------------------------------------------------------------------------


    string WaterFilterAxis::name(void) const {
      stringstream s;
      s << boost::format("WaterFilterAxis(radius=%f)") % sqrt(radius_);
      return(s.str());
    }


    vector<int> WaterFilterAxis::filter(const AtomicGroup& solv, const AtomicGroup& prot) {
      bdd_ = boundingBox(prot);

      vector<int> result(solv.size());
      AtomicGroup::const_iterator ai;
      uint j = 0;
  
      for (ai = solv.begin(); ai != solv.end(); ++ai) {
        GCoord a = (*ai)->coords();
        if (a.z() < bdd_[0][2] || a.z() > bdd_[1][2]) {
          result[j++] = 0;
          continue;
        }

        // Find the nearest point on the axis to the atom...
        a -= orig_;
        double k = (axis_ * a) / axis_.length2();
        GCoord ah = orig_ + k * axis_;
        GCoord v = (*ai)->coords() - ah;

        // Are we within the radius cutoff?
        double d = v.length2();
        if (d <= radius_)
          result[j++] = true;
        else
          result[j++] = false;
      }
      return(result);
    }

    double WaterFilterAxis::volume(void) {
      return((bdd_[1][2] - bdd_[0][2]) * PI * radius_);
    }


    vector<GCoord> WaterFilterAxis::boundingBox(const AtomicGroup& grp) {

      // Set the principal axis...
      orig_ = grp.centroid();
      vector<GCoord> axes = grp.principalAxes();
      axis_ = axes[0];
      vector<GCoord> bdd = grp.boundingBox();
  
      // Calculate the extents of the box containing the principal axis cylinder...
      double r = sqrt(radius_);
      GCoord lbd = orig_ - axis_ - GCoord(r,r,0.0);
      GCoord ubd = orig_ + axis_ + GCoord(r,r,0.0);

      // Set the z-bounds to the protein bounding box...
      lbd[2] = bdd[0][2];
      ubd[2] = bdd[1][2];

      // Replace...
      bdd[0] = lbd;
      bdd[1] = ubd;

      return(bdd);
    }

    // --------------------------------------------------------------------------------


    string WaterFilterBlob::name(void) const {
      stringstream s;
      GCoord min = blob_.minCoord();
      GCoord max = blob_.maxCoord();
      SGridpoint dim = blob_.gridDims();
      s << "WaterFilterBlob(" << dim << ":" << min << "x" << max << ")";
      return(s.str());
    }


    double WaterFilterBlob::volume(void) {
      if (vol >= 0.0)
        return(vol);

      GCoord d = blob_.gridDelta();
      double delta = d[0] * d[1] * d[2];
      long n = blob_.maxGridIndex();
      long c = 0;
      for (long i=0; i<n; ++i)
        if (blob_(i))
          ++c;

      vol = c * delta;
      return(vol);
    }

    vector<int> WaterFilterBlob::filter(const AtomicGroup& solv, const AtomicGroup& prot) {
      vector<int> result(solv.size());
      AtomicGroup::const_iterator ci;
      uint j = 0;
      for (ci = solv.begin(); ci != solv.end(); ++ci) {
        GCoord c = (*ci)->coords();
        SGridpoint probe = blob_.gridpoint(c);
        if (blob_.inRange(probe))
          result[j++] = (blob_(c) != 0);
        else
          result[j++] = 0;
      }

      return(result);
    }


    // This ignores the protein bounding box...
    vector<GCoord> WaterFilterBlob::boundingBox(const AtomicGroup& prot) {
      if (bdd_set)
        return(bdd_);
  
      SGridpoint dim = blob_.gridDims();
      SGridpoint min = dim;
      SGridpoint max(0,0,0);

      for (int k=0; k<dim[2]; ++k)
        for (int j=0; j<dim[1]; ++j)
          for (int i=0; i<dim[0]; ++i) {
            SGridpoint probe(i,j,k);
            if (blob_(probe) != 0)
              for (int x=0; x<3; ++x) {
                if (probe[x] < min[x])
                  min[x] = probe[x];
                if (probe[x] > max[x])
                  max[x] = probe[x];
              }
          }

      vector<GCoord> bdd(2);
      bdd[0] = blob_.gridToWorld(min);
      bdd[1] = blob_.gridToWorld(max);
    
      return(bdd);
    }



    // --------------------------------------------------------------------------------


    string ZClippedWaterFilter::name(void) const {
      stringstream s;

      s << boost::format("ZClippedWaterFilter(%s, %f, %f)") % WaterFilterDecorator::name() % zmin_ % zmax_;
      return(s.str());
    }


    vector<int> ZClippedWaterFilter::filter(const AtomicGroup& solv, const AtomicGroup& prot) {
      vector<int> result = WaterFilterDecorator::filter(solv, prot);

      for (uint i=0; i<result.size(); ++i)
        if (result[i]) {
          GCoord c = solv[i]->coords();
          if (c[2] < zmin_ || c[2] > zmax_)
            result[i] = 0;
        }    

      return(result);
    }


    vector<GCoord> ZClippedWaterFilter::boundingBox(const AtomicGroup& grp) {
      vector<GCoord> bdd = WaterFilterDecorator::boundingBox(grp);
      bdd[0][2] = zmin_;
      bdd[1][2] = zmax_;

      return(bdd);
    }

    // --------------------------------------------------------------------------------

    string BulkedWaterFilter::name(void) const {
      stringstream s;

      s << boost::format("BulkedWaterFilter(%s, %f, %f, %f)") % WaterFilterDecorator::name() % pad_ % zmin_ % zmax_;
      return(s.str());
    }


    vector<int> BulkedWaterFilter::filter(const AtomicGroup& solv, const AtomicGroup& prot) {
      vector<int> result = WaterFilterDecorator::filter(solv, prot);
      vector<GCoord> bdd = boundingBox(prot);

      for (uint i=0; i<result.size(); ++i)
        if (!result[i]) {
          GCoord c = solv[i]->coords();
          if ( ((c[0] >= bdd[0][0] && c[0] <= bdd[1][0]) &&
                (c[1] >= bdd[0][1] && c[1] <= bdd[1][1]) &&
                (c[2] >= bdd[0][2] && c[2] <= zmin_))
               ||
               ((c[0] >= bdd[0][0] && c[0] <= bdd[1][0]) &&
                (c[1] >= bdd[0][1] && c[1] <= bdd[1][1]) &&
                (c[2] <= bdd[1][2] && c[2] >= zmax_)) )
            result[i] = true;
        }    

      return(result);
    }


    vector<GCoord> BulkedWaterFilter::boundingBox(const AtomicGroup& grp) {
      vector<GCoord> bdd = grp.boundingBox();

      bdd[0] -= pad_;
      bdd[1] += pad_;
      return(bdd);
    }

  };


};
