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




	void AmberNetcdf::init(const char* name, const uint natoms) {
		int retval;


		retval = nc_open(name, NC_NOWRITE, &_ncid);
		if (retval)
			throw(FileOpenError(name, "", retval));    // May want to preserve code here in the future...

		// Read and validate global attributes...
		readGlobalAttributes();
		if (_conventions.empty() || _conventionVersion.empty())
			throw(FileOpenError(name, "Unable to find convention global attributes.  Is this really an Amber NetCDF trajectory?"));

		if (_conventions.find("AMBER") == std::string::npos)
			throw(FileOpenError(name, "Cannot find AMBER tag in global attribues.  Is this really an Amber NetCDF trajectory?"));

		if (_conventionVersion != std::string("1.0"))
			throw(FileOpenError(name, "Convention version is '" + _conventionVersion + "', but only 1.0 is supported for Amber NetCDF trajectories."));

		// Verify # of atoms match...
		int atom_id;
		retval = nc_inq_dimid(_ncid, "atom", &atom_id);
		if (retval)
			throw(FileOpenError(name, "Cannot read atom id"), retval);
		retval = nc_inq_dimlen(_ncid, atom_id, &_natoms);
		if (retval)
			throw(FileOpenError(name, "Cannot read atom length", retval));
		if (_natoms != natoms)
			throw(FileOpenError(name, "AmberNetcdf has different number of atoms than expected"));


		// Get nframes
		int frame_id;
		retval = nc_inq_dimid(_ncid, "frame", &frame_id);
		if (retval)
			throw(FileOpenError(name, "Cannot read frame information", retval));
		retval = nc_inq_dimlen(_ncid, frame_id, &_nframes);

		// Check for periodic cells...
		retval = nc_inq_varid(_ncid, "cell_lengths", &_cell_lengths_id);
		_periodic = !retval;

		// Any additional validation should go here...

		// Get coord-id for later use...
		retval = nc_inq_varid(_ncid, "coordinates", &_coord_id);
		if (retval)
			throw(FileOpenError(name, "Cannot get id for coordinates", retval));

		// Get velocity-id for later use...
		retval = nc_inq_varid(_ncid, "velocities", &_velocities_id);
		_velocities = !retval;


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
					throw(FileOpenError(name, "Cannot get first time point", retval));

				idx[0] = 1;
				retval = nc_get_var1_float(_ncid, time_id, idx, &t1);
				if (retval)
					throw(FileOpenError(name, "Cannot get second time point", retval));

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
	void AmberNetcdf::readRawFrame(const uint frameno)  {
		size_t start[3] = {0, 0, 0};
		size_t count[3] = {1, 1, 3};


		// Read coordinates first...
		start[0] = frameno;
		count[1] = _natoms;


		int retval = VarTypeDecider<GCoord::element_type>::read(_ncid, _coord_id, start, count, _coord_data);
		if (retval)
			throw(FileReadError(_filename, "Cannot read Amber netcdf frame (coords)", retval));

		if (_velocities)
		{
			retval = VarTypeDecider<GCoord::element_type>::read(_ncid, _velocities_id, start, count, _velocity_data);
			if (retval)
				throw(FileReadError(_filename, "Cannot read Amber netcdf frame (velocities)", retval));
		}


		// Now get box if present...
		if (_periodic) {
			start[1] = 0;
			count[1] = 3;

			retval = VarTypeDecider<GCoord::element_type>::read(_ncid, _cell_lengths_id, start, count, _box_data);
			if (retval)
				throw(FileReadError(_filename, "Cannot read Amber netcdf periodic box", retval));
		}

	}

	bool AmberNetcdf::parseFrame() {
		if (_current_frame >= _nframes)
			return(false);

		readRawFrame(_current_frame);
		return(true);
	}

	void AmberNetcdf::updateGroupCoordsImpl(AtomicGroup& g) {

		for (AtomicGroup::iterator i = g.begin(); i != g.end(); ++i) {
			uint idx = (*i)->index();
			if (idx >= _natoms)
				throw(LOOSError(_filename, **i, "Atom index into trajectory frame is out of bounds"));
			idx *= 3;
			(*i)->coords(GCoord(_coord_data[idx], _coord_data[idx+1], _coord_data[idx+2]));
		}

		if (_periodic)
			g.periodicBox(GCoord(_box_data[0], _box_data[1], _box_data[2]));
	}


	void AmberNetcdf::updateGroupVelocitiesImpl(AtomicGroup& g) {

		for (AtomicGroup::iterator i = g.begin(); i != g.end(); ++i) {
			uint idx = (*i)->index();
			if (idx >= _natoms)
				throw(LOOSError(_filename, **i, "Atom index into trajectory frame is out of bounds"));
			idx *= 3;
			(*i)->coords(GCoord(_velocity_data[idx], _velocity_data[idx+1], _velocity_data[idx+2]));
		}
	}


	std::vector<GCoord> AmberNetcdf::velocitiesImpl() const {
		std::vector<GCoord> res;
		for (uint i=0; i<_natoms; i += 3)
			res.push_back(GCoord(_velocity_data[i], _velocity_data[i+1], _velocity_data[i+2]));
		return(res);
	}


	void AmberNetcdf::readGlobalAttributes()  {

		_title = readGlobalAttribute("title");
		_application = readGlobalAttribute("application");
		_program = readGlobalAttribute("program");
		_programVersion = readGlobalAttribute("programVersion");
		_conventions = readGlobalAttribute("Conventions");
		_conventionVersion = readGlobalAttribute("ConventionVersion");
	}


	// Will return an emptry string if the attribute is not found
	// May want to consider special exception for type errors...
	std::string AmberNetcdf::readGlobalAttribute(const std::string& name) {
		size_t len;

		int retval = nc_inq_attlen(_ncid, NC_GLOBAL, name.c_str(), &len);
		if (retval)
			return(std::string());

		nc_type type;
		retval = nc_inq_atttype(_ncid, NC_GLOBAL, name.c_str(), &type);
		if (type != NC_CHAR)
			throw(FileOpenError(_filename, "Only character data is supported for global attributes", retval));

		char* buf = new char[len+1];
		retval = nc_get_att_text(_ncid, NC_GLOBAL, name.c_str(), buf);
		if (retval) {
			delete[] buf;
			throw(FileOpenError(_filename, "Cannot read attribute " + name, retval));
		}
		buf[len]='\0';

		return(std::string(buf));
	}


};
