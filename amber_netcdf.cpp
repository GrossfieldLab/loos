// (c) 2012 Tod D. Romo, Grossfield Lab, URMC

#include <amber_netcdf.hpp>
#include <AtomicGroup.hpp>

namespace loos {

  void AmberNetcdfTraj::init(const char* name, const int natoms) {
    int retval;

    retval = nc_open(name, NC_NOWRITE, &_ncid);
    if (retval)
      throw(AmberNetcdfOpenError());

    readGlobalAttributes();

    // Verify # of atoms match...
    int atom_id;
    retval = nc_inq_dimid(_ncid, "atom", &atom_id);
    if (retval)
      throw(AmberNetcdfError("Error reading atom id", retval));
    retval = nc_inq_dimlen(_ncid, atom_id, &_natoms);
    if (retval)
      throw(AmberNetcdfError("Error reading atom length", retval));
    if (_natoms != natoms)
      throw(LOOSError("AmberNetcdfTraj has different number of atoms than expected"));


    // Get nframes
    int frame_id;
    retval = nc_inq_dimid(_ncid, "frame", frame_id);
    if (retval)
      throw(AmberNetcdfError("Error reading frame information", retval));
    retval = nc_inq_dimlen(_ncid, frame_id, &_nframes);

    // Check for periodic cells...
    retval = nc_inq_varid(_ncid, "cell_lengths", &_cell_lengths_id);
    if (!retval) {
      _periodic = true;
      retval = nc_inq_vartype(_ncid, _cell_lengths_id, &_box_type);
      if (retval)
        throw(AmberNetcdfError("Error reading periodic cell data type", retval));
      if (_box_data)
        throw(LOOSError("AmberNetcdfTraj::init() internal error - box_data already exists"));
      switch(_box_type) {
      case NC_FLOAT: _box_data = new float[3]; break;
      case NC_DOUBLE: _box_data = new double[3]; break;
      default: throw(AmberNetcdfError("Only float and double supported for cell_lengths variable"));
      }
    }

    // Check angles

    // Get coord-id for later use...
    retval = nc_inq_varid(_ncid, "coordinates", &_coord_id);
    if (retval)
      throw(AmberNetcdfError("Error getting id for coordinates", retval));
    retval = nc_inq_vartype(_ncid, _coord_id, &_coord_type);
    if (retval)
      throw(AmberNetcdfError("Error getting data type for coordinates", retval));
    switch(_coord_type) {
    case NC_FLOAT: _coord_data = new float[_natoms * 3]; break;
    case NC_DOUBLE: _coord_data = new double[_natoms * 3]; break;
    default: throw(AmberNetcdfError("Only float and double supported for coordinates")); break;
    }
    
  }




  void AmberNetcdfTraj::readGlobalAttributes() {
    int retval;

    _title = readGlobalAttribute("title");
    _application = readGlobalAttribute("application");
    _program = readGlobalAttribute("program");
    _programVersion = readGlobalAttribute("programVersion");
    _conventions = readGlobalAttribute("Conventions");
    _conventionVersion = readGlobalAttribute("ConventionVersion");
  }


  // Will return an emptry string if the attribute is not found
  string AmberNetcdfTraj::readGlobalAttribute(const std::string& name) {
    size_t len;
    
    retval = nc_inq_attname(_ncid, NC_GLOBAL, name.c_str(), &len);
    if (retval)
      return(string());

    nc_type type;
    retval = nc_inq_atttype(_ncid, NC_GLOBAL, name.c_str(), &type);
    if (type != NC_CHAR)
      throw(AmberNetcdfTypeError("Only character data is supported for global attributes"));
    

    char* buf = new char[len+1];
    retval = nc_get_att_text(_ncid, NC_GLOBAL, name.c_str(), buf);
    if (retval) {
      delete[] buf;
      throw(AmberNetcdfError("Error reading attribute " + name));
    }

    return(string(buf));
  }





};
