// Exception handling...


%catches(loos::ParseError,\
	 loos::FileOpenError,\
	 loos::FileReadErrorWithLine,\
	 loos::FileReadError,\
	 loos::FileError,\
	 loos::LOOSError,\
	 std::runtime_error,\
	 std::logic_error) \
createSystem;

%catches(loos::ParseError,\
	 loos::FileOpenError,	    \
	 loos::FileReadErrorWithLine,\
	 loos::FileReadError,\
	 loos::FileError,\
	 loos::XDRDataSizeError,\
	 loos::LOOSError,\
	 std::logic_error,\
	 std::runtime_error) \
createTrajectory;


%catches(loos::NullResult,\
	 loos::ParseError,\
	 std::runtime_error,\
	 std::logic_error)\
selectAtoms;

%catches(loos::ParseError) parseRangeList;

// amber
%catches(loos::FileOpenError, loos::FileReadErrorWithLine, loos::LOOSError) Amber::Amber;
%catches(loos::FileOpenError, loos::FileReadErrorWithLine, loos::LOOSError) Amber::create;
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

// Atom
%catches(loos::LOOSError) Atom::checkUserBits;
%catches(loos::LOOSError) Atom::setProperty;
%catches(loos::LOOSError) Atom::clearProperty;
%catches(loos::UnsetProperty) Atom::charge;
%catches(loos::LOOSError) Atom::deleteBond;
%catches(loos::UnsetProperty) Atom::getBonds;



// AtomicGroup
%catches(std::out_of_range) AtomicGroup::rangeCheck;
%catches(std::out_of_range) AtomicGroup::getAtom;
%catches(std::out_of_range) AtomicGroup::operator[];
%catches(loos::LOOSError) AtomicGroup::deleteAtom;
%catches(std::out_of_range) AtomicGroup::subset;
%catches(std::out_of_range) AtomicGroup::excise;
%catches(loos::LOOSError) AtomicGroup::reimage;
%catches(loos::LOOSError) AtomicGroup::reimageByAtom;
%catches(loos::LOOSError) AtomicGroup::mergeImage;
%catches(loos::LOOSError) AtomicGroup::atomOrderMapFrom;
%catches(loos::LOOSError) AtomicGroup::copyMappedCoordinatesFrom;

%catches(loos::LOOSError) AtomicGroup::rmsd;

%catches(loos::NumericalError) AtomicGroup::momentsOfInertia;
%catches(loos::NumericalError) AtomicGroup::principalAxes;
%catches(loos::NumericalError) AtomicGroup::superposition;
%catches(loos::NumericalError) AtomicGroup::alignOnto;


// ccpdb
%catches(loos::FileReadError, loos::LOOSError) CCPDB::parseFrame;
%catches(loos::FileOpenError, loos::FileReadError, loos::LOOSError) CCPDB::CCPDB;
%catches(loos::FileOpenError, loos::LOOSError) CCPDB::create;

// charmm
%catches(loos::FileOpenError, loos::FileReadError) CHARMM::CHARMM;
%catches(loos::FileReadError) CHARMM::read;

// dcd
%catches(loos::FileOpenError, loos::FileReadError, loos::FileError, loos::LOOSError) DCD::DCD;
%catches(loos::FileOpenError, loos::FileReadError, loos::FileError, loos::LOOSError) DCD::create;
%catches(loos::FileReadError, loos::LOOSError) DCD::parseFrame;

// dcdwriter
%catches(std::logic_error) DCDWriter::setHeader;
%catches(std::logic_error) DCDWriter::setTitles;
%catches(std::logic_error) DCDWriter::setTitle;
%catches(std::logic_error) DCDWriter::addTitle;
%catches(std::logic_error) DCDWriter::setComments;
%catches(loos::FileWriteError, loos::LOOSError) DCDWriter::writeFrame;
%catches(loos::FileOpenError, loos::FileReadError, loos::FileError, loos::LOOSError) DCDWriter::prepareToAppend;

// gro

%catches(loos::FileOpenError, loos::FileReadError, loos::ParseError, loos::LOOSError) Gromacs::Gromacs;
%catches(loos::FileOpenError, loos::FileReadError, loos::ParseError, loos::LOOSError) Gromacs::create;

// MDTraj 
%catches(loos::FileOpenError, loos::FileReadError, loos::LOOSError) MDTraj::MDTraj;
// TODO: there are some H5 exceptions I should probably catch here
%catches(loos::FileReadError, loos::LOOSError, H5::FileIException) MDTraj::read;


// pdb_remarks

%catches(std::range_error) Remarks::get;
%catches(std::range_error) Remarks::erase;
%catches(std::range_error) Remarks::rangeCheck;


// pdb

%catches(loos::ParseError) PDB::parseAtomRecord;
%catches(loos::LOOSError, std::logic_error) PDB::atomAsString;
%catches(loos::LOOSError) PDB::findAtom;
%catches(loos::LOOSError, loos::ParseError) PDB::parseConectRecord;
%catches(loos::ParseError) PDB::parseCryst1Record;
%catches(loos::FileReadError) PDB::read;
%catches(loos::LOOSError) PDB::FormatConectRecords;
%catches(loos::FileOpenError, loos::FileReadError, loos::LOOSError) PDB::PDB;
%catches(loos::FileOpenError, loos::FileReadError, loos::LOOSError) PDB::create;


// pdbtraj
%catches(loos::FileReadError, loos::LOOSError) PDBTraj::parseFrame;
%catches(loos::FileOpenError) PDBTraj::PDBTraj;

// psf
%catches(loos::FileOpenError, loos::FileReadError, loos::FileError, loos::LOOSError) PSF::PSF;
%catches(loos::FileOpenError, loos::FileReadError, loos::FileError, loos::LOOSError) PSF::create;
%catches(loos::FileReadError, loos::FileError, loos::LOOSError) PSF::read;

// TinkerArc
%catches(loos::FileOpenError, loos::FileReadError, loos::FileError, loos::LOOSError) TinkerArc::TinkerArc;
%catches(loos::FileOpenError, loos::FileReadError, loos::FileError, loos::LOOSError) TinkerArc::create;
%catches(loos::FileReadError, loos::LOOSError) TinkerArc::parseFrame;

// tinkerxyz
%catches(loos::FileOpenError, loos::FileReadError, loos::ParseError, loos::LOOSError) TinkerXYZ::TinkerXYZ;
%catches(loos::FileOpenError, loos::FileReadError, loos::ParseError, loos::LOOSError) TinkerXYZ::create;
%catches(loos::FileReadError, loos::FileError, loos::LOOSError) TinkerXYZ::read;

// trajectory
%catches(loos::ParseError,\
	 loos::FileOpenError,	    \
	 loos::FileReadErrorWithLine,\
	 loos::FileReadError,\
	 loos::FileError,\
	 loos::XDRDataSizeError,\
	 loos::LOOSError,\
	 std::logic_error,\
	 std::runtime_error) \
Trajectory::Trajectory;

%catches(loos::ParseError,\
	 loos::FileOpenError,	    \
	 loos::FileReadErrorWithLine,\
	 loos::FileReadError,\
	 loos::FileError,\
	 loos::XDRDataSizeError,\
	 loos::LOOSError,\
	 std::logic_error,\
	 std::runtime_error) \
Trajectory::readFrame;

%catches(loos::ParseError,\
	 loos::FileOpenError,	    \
	 loos::FileReadErrorWithLine,\
	 loos::FileReadError,\
	 loos::FileError,\
	 loos::XDRDataSizeError,\
	 loos::LOOSError,\
	 std::logic_error,\
	 std::runtime_error) \
Trajectory::parseFrame;

%catches(loos::FileError, std::range_error) Trajectory::seekNextFrame;
%catches(loos::FileError, std::range_error) Trajectory::seekFrame;
%catches(loos::FileError) Trajectory::rewind;
%catches(loos::LOOSError, std::range_error) Trajectory::updateGroupCoords;



// trajwriter
%catches(loos::FileOpenError, loos::FileReadError, loos::FileError, loos::LOOSError) TrajectoryWriter::TrajectoryWriter;
%catches(loos::FileWriteError, loos::FileError, loos::LOOSError) TrajectoryWriter::writeFrame;
%catches(std::logic_error) TrajectoryWriter::setComments;

// utils.cpp

%catches(loos::ParseError, loos::LOOSError) selectAtoms;
%catches(std::logic_error) parseStringAsHybrid36;
%catches(std::logic_error, loos::LOOSError) hybrid36AsString;

%catches(std::ios_base::failure) getNextLine;
%catches(loos::LOOSError) invocationHeader;


// xtc
%catches(loos::FileOpenError, loos::FileReadError, loos::FileError, loos::XDRDataSizeError, loos::LOOSError) XTC::XTC;
%catches(loos::FileOpenError, loos::FileReadError, loos::FileError, loos::XDRDataSizeError, loos::LOOSError) XTC::create;

// xtcwrite

%catches(loos::FileReadError, loos::FileOpenError, loos::FileError, loos::LOOSError) XTCWriter::XTCWriter;
%catches(loos::FileWriteError, loos::LOOSError) XTCWriter::writeFame;
