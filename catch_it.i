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


// amber
%catches(loos::FileOpenError, loos::FileReadErrorWithLine) Amber::Amber;
%catches(loos::FileReadErrorWithLine) Amber::parseFormat;
%catches(loos::FileReadErrorWithLine) Amber::parseCharges;
%catches(loos::FileReadErrorWithLine) Amber::parseMasses;
%catches(loos::FileReadErrorWithLine) Amber::parseResidueLabels;
%catches(loos::FileReadErrorWithLine) Amber::parseResiduePointers;
%catches(std::runtime_error) Amber::assignResidues;
%catches(loos::FileReadErrorWithLine) Amber::parseBonds;
%catches(std::logic_error) Amber::parsePointers;
%catches(loos::FileReadErrorWithLine) Amber::parseTitle;
%catches(loos::FileReadErrorWithLine) Amber::parseAtomNames;
%catches(loos::FileReadErrorWithLine) Amber::parseAmoebaRegularBondNumList;
%catches(loos::FileReadErrorWithLine) Amber::parseAmoebaRegularBondList;
%catches(loos::FileReadErrorWithLine, std::logic_error) Amber::read;


// amber_netcdf
%catches(loos::FileOpenError) AmberNetcdf::init;
%catches(loos::FileReadError) AmberNetcdf::readRawFrame;
%catches(loos::FileReadError) AmberNetcdf::parseFrame();
%catches(loos::LOOSError) AmberNetcdf::updateGroupCoordsImpl;
%catches(loos::FileOpenError) AmberNetcdf::readGlobalAttribute;
%catches(loos::FileOpenError) Ambernetcdf::readGlobalAttributes;
%catches(loos::FileError) Ambernetcdf::seekFrameImpl;

// amber_rst
%catches(loos::FileOpenError) AmberRst::AmberRst;
%catches(loos::FileReadError) AmberRst::parseFrame;
%catches(loos::LOOSError) AmberRst::updateGroupCoordsImpl;

// amber_traj
%catches(loos::FileOpenError) AmberTraj::AmberTraj;
%catches(loos::FileError) AmberTraj::parseFrame;
%catches(loos::FileError) AmberTraj::seekFrameImpl;
%catches(loos::LOOSError) AmberTraj::updateGroupCoordsImpl;


// utils.cpp

%catches(loos::ParseError) selectAtoms;
%catches(std::logic_error) parseStringAsHybrid36;
%catches(std::logic_error, loos::LOOSError) hybrid36AsString;
