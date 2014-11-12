// Exception handling...


%catches(loos::FileParseError,\
	 loos::BadConnectivityError,\
	 loos::FileOpenError,\
	 loos::FileReadError,\
	 loos::LOOSError,     \
	 std::runtime_error,\
	 std::logic_error) \
createSystem;

%catches(loos::AmberNetcdfOpenError,\
	 loos::AmberNetcdfError,\
	 loos::FileOpenError,\
	 loos::FileReadError,\
	 loos::LOOSError,\
	 loos::TrajectoryReadError,\
	 std::logic_error,\
	 std::runtime_error) \
createTrajectory;


%catches(loos::NullResult,\
	 loos::ParseError,\
	 std::runtime_error,\
	 std::logic_error)\
selectAtoms;


// amber_netcdf
%catches(loos::FileOpenError) AmberNetcdf::init;
%catches(loos::FileReadError) AmberNetcdf::readRawFrame;
%catches(loos::FileReadError) AmberNetcdf::parseFrame();
%catches(loos::LOOSError) AmberNetcdf::updateGroupCoordsImpl;
%catches(loos::FileOpenError) AmberNetcdf::readGlobalAttribute;
%catches(loos::FileOpenError) Ambernetcdf::readGlobalAttributes;

// utils.cpp

%catches(loos::ParseError) selectAtoms;
%catches(std::logic_error) parseStringAsHybrid36;
%catches(std::logic_error, loos::LOOSError) hybrid36AsString;
