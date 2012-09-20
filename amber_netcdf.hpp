// (c) 2012 Tod D. Romo, Grossfield Lab, URMC

#if !defined(LOOS_AMBER_NETCDF_HPP)
#define LOOS_AMBER_TRAJ_HPP


#include <iostream>
#include <string>
#include <netcdf.h>

#include <loos_defs.hpp>
#include <Coord.hpp>
#include <Trajectory.hpp>
#include <exceptions.hpp>

namespace loos {

  struct AmberNetcdfError : public LOOSError {
    explicit AmberNetcdfError(const std::string& msg) : LOOSError(msg) { }
    explicit AmberNetcdfError(const std::string& msg, const int retval) {
      std::stringstream ss;
      ss << msg << " with error #" << retval;
      _msg = ss.str();
    }

  };

  struct AmberNetcdfOpenError : public AmberNetcdfError {
    explicit AmberNetcdfOpenError() : AmberNetcdfError("Error opening Amber NetCDF file") { }
  };

  struct AmberNetcdfTypeError : public AmberNetcdfError {
    explicit AmberNetcdfTypeError(const std::string msg) : AmberNetcdfError(msg) { }
  };


  
  class AmberNetcdf : public Trajectory {
  public:

    
    explicit AmberNetcdf(const std::string& s, const int na) :
      _coord_data(0),
      _box_data(0),
      _periodic(false),
      _current_frame(0)
    {
      cached_first = false;
      init(s.c_str(), na);
    }

    explicit AmberNetcdf(const char* p, const int na) :
      _coord_data(0),
      _box_data(0),
      _periodic(false),
      _current_frame(0)
    {
      cached_first = false;
      init(p, na);
    }

    ~AmberNetcdf() {
      int retval = nc_close(_ncid);
      if (retval)
        throw(AmberNetcdfError("Error while closing netcdf file", retval));

      if (_coord_data) {
        switch(_coord_type) {
        case NC_FLOAT: delete[] static_cast<float*>(_coord_data); break;
        case NC_DOUBLE: delete[] static_cast<double*>(_coord_data); break;
        }
      }

      if (_box_data) {
        switch(_box_type) {
        case NC_FLOAT: delete[] static_cast<float*>(_box_data); break;
        case NC_DOUBLE: delete[] static_cast<double*>(_box_data); break;
        }
      }
    }

    uint natoms() const { return(_natoms); }
    uint nframes() const { return(_nframes); }
    float timestep() const { return(0.0); }   // Fix!
    
    bool hasPeriodicBox() const { return(_periodic); }
    GCoord periodicBox() const { return(GCoord(_box_data[0], _box_data[1], _box_data[2])); }

    std::vector<GCoord> coords() {
      std::vector<GCoord> res;
      for (uint i=0; i<_natoms; i += 3)
        res.push_back(GCoord(_coord_data[i], _coord_data[i+1], _coord_data[i+2]));
      return(res);
    }

  private:
    void init(const char* name, const int natoms);


  private:
    void* _coord_data;
    void* _box_data;
    bool _periodic;
    uint _current_frame;
    int _ncid;
    size_t _nframes;
    size_t _natoms;
    int _coord_id;
    nc_type _coord_type, _box_type;
    size_t _coord_size;
    int _cell_lengths_id;
    std::string _title, _application, _program, _programVersion, _conventions, _conventionVersion;
    std::vector<GCoord> _frame;
    GCoord _box;
  };


}



#endif
