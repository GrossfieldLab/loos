// (c) 2012 Tod D. Romo, Grossfield Lab, URMC

#include <amber_netcdf.hpp>
#include <AtomicGroup.hpp>

namespace loos {


  bool isFileNetCDF(const std::string& fname) {
    std::ifstream ifs(fname.c_str());

    char buf[4];
    ifs.read(buf, 4);
    return (buf[0] == 'C' && buf[1] == 'D' && buf[2] == 'F' && (buf[3] = 0x01 || buf[3] == 0x02));
  }




  void AmberNetcdf::init(const char* name, const uint natoms) throw(AmberNetcdfOpenError, AmberNetcdfError){
    int retval;

    retval = nc_open(name, NC_NOWRITE, &_ncid);
    if (retval)
      throw(AmberNetcdfOpenError());

    // Read and validate global attributes...
    readGlobalAttributes();
    if (_conventions.empty() || _conventionVersion.empty())
      throw(AmberNetcdfError("Unable to find convention global attributes.  Is this really an Amber NetCDF trajectory?"));

    if (_conventions.find("AMBER") == std::string::npos)
      throw(AmberNetcdfError("Cannot find AMBER tag in global attribues.  Is this really an Amber NetCDF trajectory?"));
      
    if (_conventionVersion != std::string("1.0"))
      throw(AmberNetcdfError("Convention version is '" + _conventionVersion + "', but only 1.0 is supported for Amber NetCDF trajectories."));

    // Verify # of atoms match...
    int atom_id;
    retval = nc_inq_dimid(_ncid, "atom", &atom_id);
    if (retval)
      throw(AmberNetcdfError("Error reading atom id", retval));
    retval = nc_inq_dimlen(_ncid, atom_id, &_natoms);
    if (retval)
      throw(AmberNetcdfError("Error reading atom length", retval));
    if (_natoms != natoms)
      throw(AmberNetcdfError("AmberNetcdf has different number of atoms than expected"));


    // Get nframes
    int frame_id;
    retval = nc_inq_dimid(_ncid, "frame", &frame_id);
    if (retval)
      throw(AmberNetcdfError("Error reading frame information", retval));
    retval = nc_inq_dimlen(_ncid, frame_id, &_nframes);

    // Check for periodic cells...
    retval = nc_inq_varid(_ncid, "cell_lengths", &_cell_lengths_id);
    _periodic = !retval;

    // Any additional validation should go here...

    // Get coord-id for later use...
    retval = nc_inq_varid(_ncid, "coordinates", &_coord_id);
    if (retval)
      throw(AmberNetcdfError("Error getting id for coordinates", retval));

    // Attempt to determine timestep by looking at dT between frames 1 & 2
    if (_nframes >= 2) {
      int time_id;
      retval = nc_inq_varid(_ncid, "time", &time_id);
      if (!retval) {
        float t0, t1;
        size_t idx[1];

        idx[0] = 0;
        retval = nc_get_var1_float(_ncid, time_id, idx, &t0);
        if (retval)
          throw(AmberNetcdfError("Error getting first time point", retval));

        idx[0] = 1;
        retval = nc_get_var1_float(_ncid, time_id, idx, &t1);
        if (retval)
          throw(AmberNetcdfError("Error getting second time point", retval));

        // Assume units are in picoseconds
        _timestep = (t1-t0)*1e-12;
      }
    }


    // Now cache the first frame...
    readRawFrame(0);
    cached_first = true;
    
  }


  // Given a frame number, read the coord data into the internal array
  // and retrieve the corresponding periodic box (if present)
  void AmberNetcdf::readRawFrame(const uint frameno) throw(AmberNetcdfError) {
    size_t start[3] = {0, 0, 0};
    size_t count[3] = {1, 1, 3};


    // Read coordinates first...
    start[0] = frameno;
    count[1] = _natoms;


    int retval = VarTypeDecider<GCoord::element_type>::read(_ncid, _coord_id, start, count, _coord_data);
    if (retval)
      throw(AmberNetcdfError("Error while reading Amber netcdf frame", retval));

    // Now get box if present...
    if (_periodic) {
      start[1] = 0;
      count[1] = 3;

      retval = VarTypeDecider<GCoord::element_type>::read(_ncid, _cell_lengths_id, start, count, _box_data);
      if (retval)
        throw(AmberNetcdfError("Error while reading Amber netcdf periodic box", retval));
      
    }
    
  }


  void AmberNetcdf::seekNextFrameImpl() {
    ++_current_frame;
  }

  void AmberNetcdf::seekFrameImpl(const uint i) {
    _current_frame = i;
  }

  void AmberNetcdf::rewindImpl() {
    _current_frame = 0;
  }

  bool AmberNetcdf::parseFrame() throw(AmberNetcdfError) {
    if (_current_frame >= _nframes)
      return(false);

    readRawFrame(_current_frame);
    return(true);
  }

  void AmberNetcdf::updateGroupCoordsImpl(AtomicGroup& g) throw(LOOSError) {

    for (AtomicGroup::iterator i = g.begin(); i != g.end(); ++i) {
      uint idx = (*i)->index();
      if (idx >= _natoms)
        throw(LOOSError(**i, "Atom index into trajectory frame is out of bounds"));
      idx *= 3;
      (*i)->coords(GCoord(_coord_data[idx], _coord_data[idx+1], _coord_data[idx+2]));
    }
    
    if (_periodic)
      g.periodicBox(GCoord(_box_data[0], _box_data[1], _box_data[2]));
  }


  void AmberNetcdf::readGlobalAttributes() throw(AmberNetcdfError) {

    _title = readGlobalAttribute("title");
    _application = readGlobalAttribute("application");
    _program = readGlobalAttribute("program");
    _programVersion = readGlobalAttribute("programVersion");
    _conventions = readGlobalAttribute("Conventions");
    _conventionVersion = readGlobalAttribute("ConventionVersion");
  }


  // Will return an emptry string if the attribute is not found
  std::string AmberNetcdf::readGlobalAttribute(const std::string& name) throw(AmberNetcdfError) {
    size_t len;
    
    int retval = nc_inq_attlen(_ncid, NC_GLOBAL, name.c_str(), &len);
    if (retval)
      return(std::string());

    nc_type type;
    retval = nc_inq_atttype(_ncid, NC_GLOBAL, name.c_str(), &type);
    if (type != NC_CHAR)
      throw(AmberNetcdfTypeError("Only character data is supported for global attributes"));
    

    char* buf = new char[len+1];
    retval = nc_get_att_text(_ncid, NC_GLOBAL, name.c_str(), buf);
    buf[len]='\0';
    if (retval) {
      delete[] buf;
      throw(AmberNetcdfError("Error reading attribute " + name));
    }

    return(std::string(buf));
  }





};
