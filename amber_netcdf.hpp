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


    struct SetEdgesForCoords {
      static const uint ndims = 3;

      SetEdgesForCoords(const uint natoms) : _natoms(natoms) { }

      void setStart(size_t s[], const uint frame) const {
        s[0] = frame;
        s[1] = 0;
        s[2] = 0;
      }

      void setCount(size_t s[]) const {
        s[0] = 1;
        s[1] = _natoms;
        s[2] = 3;
      }

      size_t size() const { return(_natoms * 3); }


      uint _natoms;
    };


    struct SetEdgesForBoxes {
      static const uint ndims = 2;

      void setStart(size_t s[], const uint frame) const {
        s[0] = frame;
        s[1] = 0;
      }

      void setCount(size_t s[]) const {
        s[0] = 1;
        s[1] = 3;
      }

      size_t size() const { return(3); }

    };



    template<class SETTER>
    struct VariableWrapper {
      VariableWrapper(const int ncid, const int varid, const SETTER& setter) :
        _ncid(ncid),
        _varid(varid),
        _setter(setter),
        _data(0)

      {
        init();
      }

      VariableWrapper(const SETTER& setter) :
        _ncid(-1),
        _varid(-1),
        _setter(setter),
        _data(0)
      { }

      void ncid(const int i) { _ncid = i; }
      int ncid() const { return(_ncid); }
      
      void varid(const int i) { _varid = i; }
      int varid() const { return(_varid); }


      void freeSpace() {
        if (data != 0)
          switch(_type) {
          case NC_BYTE: delete[] (static_cast<unsigned char*>(_data)); break;
          case NC_CHAR: delete[] (static_cast<char*>(_data)); break;
          case NC_SHORT: delete[] (static_cast<short*>(_data)); break;
          case NC_INT: delete[] (static_cast<int*>(_data)); break;
          case NC_FLOAT: delete[] (static_cast<float*>(_data)); break;
          case NC_DOUBLE: delete[] (static_cast<double*>(_data)); break;
          }
      }


      void init() {
        if (data)
          freeSpace();

        // Get type...
        int retval = nc_inq_vartype(_ncid, _varid, &_type);

        // Allocate space
        switch(_type) {
        case NC_BYTE: _data = new unsigned char[_setter.size()]; break;
        case NC_CHAR: _data = new char[_setter.size()]; break;
        case NC_SHORT: _data = new short[_setter.size()]; break;
        case NC_INT: _data=new int[_setter.size()]; break;
        case NC_FLOAT: _data=new float[_setter.size()]; break;
        case NC_DOUBLE: _data=new double[_setter.size()]; break;
        }
      }


      ~VariableWrapper() {
        freeSpace();
      }


      template<typename T>
      T getValueAt(const uint i) {
        T x;

        switch(_type) {
        case NC_BYTE: x = (static_cast<unsigned char*>(_data))[i]; break;
        case NC_CHAR: x = (static_cast<char*>(_data))[i]; break;
        case NC_SHORT: x = (static_cast<short*>(_data))[i]; break;
        case NC_INT: x = (static_cast<int*>(_data))[i]; break;
        case NC_FLOAT: x = (static_cast<float*>(_data))[i]; break;
        case NC_DOUBLE: x = (static_cast<double*>(_data))[i]; break;
        }
        return(x);
      }
                        

      nc_type type() const { return(_type); }


      bool readFrame(const uint frame) {
        _setter.setStart(_start, frame);
        _setter.setCount(_count);

        int retval = nc_get_vara(_ncid, _varid, _start, _count, _data);
        if (retval)
          cerr << "Internal error - " << retval << endl;
        return(retval == 0);
      }


      int _ncid, _varid;
      const SETTER _setter;
      void* _data;
      nc_type _type;
      size_t _start[SETTER::ndims];
      size_t _count[SETTER::ndims];
    };



  public:

    
    explicit AmberNetcdf(const std::string& s, const int na) :
      _coord_wrapper(SetEdgesForCoords(na)),
      _box_wrapper(SetEdgesForBoxes()),
      _periodic(false),
      _current_frame(0)
    {
      cached_first = false;
      init(s.c_str(), na);
    }

    explicit AmberNetcdf(const char* p, const int na) :
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
    VariableWrapper<SetEdgesForCoords> _coord_wrapper;
    VariableWrapper<SetEdgesForBoxes> _box_wrapper;
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
