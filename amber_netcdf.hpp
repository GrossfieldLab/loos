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


  
  class AmberNetcdfTraj : public Trajectory {
  public:

    
    explicit AmberNetcdfTraj(const std::string& s, const int na) :
      Trajectory(s),
      _coord_data(0),
      _box_data(0),
      _periodic(false)
    {
      ifs()->close(); // Subvert normal mechanism...
      init(s.c_str(), na);
    }

    explicit AmberNetcdfTraj(const char* p, const int na) :
      Trajectory(s),
      _coord_data(0),
      _box_data(0),
      _periodic(false)
    {
      ifs()->close();
      init(p, na);
    }

  private:
    void init(const char* name, const int natoms);


  private:
    void* _coord_data;
    void* _box_data;
    bool _periodic;
    int _ncid;
    size_t _nframes;
    size_t _natoms;
    int _coord_id;
    nc_type _coord_type, _box_type;
    int _cell_lengths_id;
    std::string _title, _application, _program, _programVersion, _conventions, _conventionVersion;
    std::vector<GCoord> _frame;
    GCoord _box;
  };


}
