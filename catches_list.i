// Exception handling...


%catches(loos::FileParseError,\
	 loos::LOOSError,     \
	 std::runtime_error,\
	 std::logic_error) \
createSystem;


