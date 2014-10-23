// Exception handling...


%catches(loos::FileParseError,\
	 loos::BadConnectivityError,\
	 loos::LOOSError,     \
	 std::runtime_error,\
	 std::logic_error) \
createSystem;

%catches(loos::AmberNetcdfOpenError,\
	 loos::AmberNetcdfError,\
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



