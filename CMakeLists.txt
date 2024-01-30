cmake_minimum_required(VERSION 3.18)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/modules/" ${CMAKE_MODULE_PATH})
cmake_policy(SET CMP0115 NEW)

include(CheckTypeSize)


project(
  LOOS
  VERSION 4.1.0
  DESCRIPTION "The Lightweight Object Oriented Structural analysis library/toolkit"
  )

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_PYLOOS "Build Python interface and install python components" ON)

set(CMAKE_CXX_STANDARD          14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(MACOSX TRUE)
endif()

check_type_size("ulong" SIZEOF_ULONG)
if(NOT ${HAVE_SIZEOF_ULONG})
  message("no ulong found")
  add_compile_definitions(REQUIRES_ULONG)
endif()

find_package(Boost REQUIRED COMPONENTS
  regex program_options json)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()


find_package(NetCDF REQUIRED)
include_directories(${NetCDF_INCLUDE_DIRS})

find_package(gemmi)

set(BLA_VENDOR OpenBLAS)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

find_package(HDF5 REQUIRED COMPONENTS CXX C HL)

if (HDF5_FOUND)
  include_directories(${HDF5_INCLUDE_DIRS})
endif()


find_package(SWIG 4.0 COMPONENTS python)
if(SWIG_FOUND)
  include(${SWIG_USE_FILE})
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)
else()
  message(WARNING "SWIG Python bindings will not be built")
endif()

include(GNUInstallDirs)

set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
if(NOT APPLE)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_subdirectory(src)
add_subdirectory(Tools)
add_subdirectory(Packages)
add_subdirectory(share)

# Doxygen handling...
find_package(Doxygen COMPONENTS dot)
if(DOXYGEN_FOUND)
  message(STATUS "Documentation can be built using --target=docs")
  set(DOXYGEN_EXCLUDE_PATTERNS "*.cc;*.hh")
  doxygen_add_docs(
    docs
    ${PROJECT_SOURCE_DIR}
    )
else()
  message(WARNING "Documentation cannot be built...requires doxygen and graphviz")
endif()

configure_file(loos_setup.sh.in loos_setup.sh)
configure_file(loos_setup.csh.in loos_setup.csh)
install(PROGRAMS ${CMAKE_BINARY_DIR}/loos_setup.sh ${CMAKE_BINARY_DIR}/loos_setup.csh DESTINATION bin)
#install(TARGETS loos pyloos)
