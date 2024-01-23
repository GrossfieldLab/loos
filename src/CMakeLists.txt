set(loos_private_files
  grammar.cc
  grammar.hh
  location.hh
  position.hh
  scanner.cc
  stack.hh
)
set(loos_hdr_files
  Atom.hpp
  AtomicGroup.hpp
  AtomicNumberDeducer.hpp
  Coord.hpp
  Fmt.hpp
  FormFactor.hpp
  FormFactorSet.hpp
  Geometry.hpp
  HBondDetector.hpp
  Kernel.hpp
  KernelActions.hpp
  KernelStack.hpp
  KernelValue.hpp
  LineReader.hpp
  LoosLexer.hpp
  Matrix.hpp
  Matrix44.hpp
  MatrixIO.hpp
  MatrixImpl.hpp
  MatrixOps.hpp
  MatrixOrder.hpp
  MatrixRead.hpp
  MatrixStorage.hpp
  MatrixUtils.hpp
  MatrixWrite.hpp
  MultiTraj.hpp
  OptionsFramework.hpp
  Parser.hpp
  ParserDriver.hpp
  PeriodicBox.hpp
  ProgressCounters.hpp
  ProgressTriggers.hpp
  RnaSuite.hpp
  Selectors.hpp
  Simplex.hpp
  StreamWrapper.hpp
  TimeSeries.hpp
  Trajectory.hpp
  UniformWeight.hpp
  UniqueStrings.hpp
  Weights.hpp
  WeightsFromFile.hpp
  XForm.hpp
  alignment.hpp
  amber.hpp
  amber_netcdf.hpp
  amber_rst.hpp
  amber_traj.hpp
  ccpdb.hpp
  charmm.hpp
  cryst.hpp
  dcd.hpp
  dcd_utils.hpp
  dcdwriter.hpp
  ensembles.hpp
  exceptions.hpp
  gro.hpp
  index_range_parser.hpp
  loos.hpp
  loos_defs.hpp
  loos_timer.hpp
  mdtraj.hpp
  mdtrajtraj.hpp
  mmcif.hpp
  pdb.hpp
  pdb_remarks.hpp
  pdbtraj.hpp
  psf.hpp
  sfactories.hpp
  sorting.hpp
  tinker_arc.hpp
  tinkerxyz.hpp
  trajwriter.hpp
  trr.hpp
  utils.hpp
  utils_random.hpp
  utils_structural.hpp
  version.hpp
  xdr.hpp
  xtc.hpp
  xtcwriter.hpp
)

set(loos_src_files
  AG_linalg.cpp
  AG_numerical.cpp
  Atom.cpp
  AtomicGroup.cpp
  AtomicNumberDeducer.cpp
  Fmt.cpp
  FormFactor.cpp
  FormFactorSet.cpp
  Geometry.cpp
  HBondDetector.cpp
  Kernel.cpp
  KernelActions.cpp
  KernelStack.cpp
  KernelValue.cpp
  LineReader.cpp
  MatrixOps.cpp
  MultiTraj.cpp
  OptionsFramework.cpp
  ProgressCounters.cpp
  ProgressTriggers.cpp
  RnaSuite.cpp
  Selectors.cpp
  UniformWeight.cpp
  Weights.cpp
  WeightsFromFile.cpp
  XForm.cpp
  alignment.cpp
  amber.cpp
  amber_netcdf.cpp
  amber_rst.cpp
  amber_traj.cpp
  ccpdb.cpp
  charmm.cpp
  dcd.cpp
  dcdwriter.cpp
  ensembles.cpp
  gro.cpp
  index_range_parser.cpp
  mdtraj.cpp
  mdtrajtraj.cpp
  mmcif.cpp
  pdb.cpp
  pdb_remarks.cpp
  pdbtraj.cpp
  psf.cpp
  sfactories.cpp
  tinker_arc.cpp
  tinkerxyz.cpp
  trr.cpp
  utils.cpp
  utils_random.cpp
  utils_structural.cpp
  version.cpp
  xtc.cpp
  xtcwriter.cpp
)

string(TIMESTAMP LOOS_BUILD_DATE "%Y%m%d")
string(TIMESTAMP LOOS_BUILD_TIMESTAMP "%Y%m%d %H:%M:%S")
configure_file(version.hpp.in version.hpp)
configure_file(version.cpp.in version.cpp)

add_library(loos ${loos_src_files} ${loos_hdr_files} ${loos_private_files})

target_include_directories(loos PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(loos PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
#target_compile_features(loos PUBLIC cxx_std_14)

# Uncomment to install header files
#set_target_properties(loos PROPERTIES PUBLIC_HEADER "${loos_hdr_files}")

target_link_libraries(loos Boost::program_options
  Boost::regex
  Boost::json
  NetCDF::NetCDF
  BLAS::BLAS
  LAPACK::LAPACK
  ${HDF5_CXX_LIBRARIES}
  gemmi::gemmi_cpp)


if(SWIG_FOUND AND BUILD_PYLOOS)
  set(CMAKE_SWIG_FLAGS -doxygen)
  set_property(SOURCE loos.i PROPERTY CPLUSPLUS ON)
  set_source_files_properties("loos.i" PROPERTIES USE_SWIG_DEPENDENCIES TRUE)   # Needed for Ninja

  swig_add_library(pyloos TYPE SHARED LANGUAGE python
    SOURCES
    loos.i
    )

  add_custom_target(loos_python ALL
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/loos ${CMAKE_CURRENT_BINARY_DIR}/pyloos
  )

  add_dependencies(loos_python pyloos)

  if(MACOSX)
    set_target_properties(pyloos PROPERTIES SUFFIX ".so")
  endif()

  target_link_libraries(pyloos Python3::Python Python3::NumPy loos)
  if(${MACOSX})
    target_link_options(pyloos PRIVATE -flat_namespace -undefined suppress)
  endif()

  add_custom_command(TARGET pyloos POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/_pyloos.so ${CMAKE_CURRENT_BINARY_DIR}/pyloos/src/loos/_pyloos.so
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/loos.py ${CMAKE_CURRENT_BINARY_DIR}/pyloos/src/loos/loos.py
 )

 install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install ${CMAKE_CURRENT_BINARY_DIR}/pyloos)")

 install(TARGETS pyloos DESTINATION lib)


endif()

install(TARGETS loos DESTINATION lib)
