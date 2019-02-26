# Operating System Compatibility

                   | LOOS    | PyLOOS  |
Operating System   | Support | Support | Notes
----------------   | ------- | ------- | -----
Fedora 17          | yes     | yes     | Unsupported
Fedora 18          | yes     | yes     | Deprecated
Fedora 19          | yes     | yes     |
Fedora 20          | yes     | yes     |
Fedora 21          | yes     | yes     |
Fedora 22          | yes     | yes     |
Fedora 23          | yes     | yes     |
Fedora 24          | yes     | yes     |
Fedora 25          | yes     | yes     |
Ubuntu 12.04 LTS   | yes     | yes     |
Ubuntu 14.04 LTS   | yes     | yes     |
Ubuntu 15.04       | yes     | yes     |
Ubuntu 15.10       | yes     | yes     |
Ubuntu 16.04 LTS   | yes     | yes     |
Debian 7.8         | yes     | yes     |
Debian 8.1         | yes     | yes     |
Centos 6.7         | yes     | yes     | See OS notes
Centos 7           | yes     | yes     | See OS notes
OpenSUSE 12        | yes     | yes     | See OS notes
OpenSUSE 13        | yes     | yes     |
Manjaro 0.8        | yes     | yes     | Unsupported
Mint 17            | yes     | yes     | Unsupported
Mint 18            | yes     | yes     |
Slackware 14.1     | yes     | yes     | Unsupported
Windows 7 (Cygwin) | yes     | no      | Unsupported
MacOS 10.11        | yes     | yes     | See OS notes


* Deprecated: We used to support this configuration, but no longer test it.  It may still work.
* Unsupported: We have built LOOS in the past using this configuration, but do not
  regularly test it and provide no direct support for using it.

As of LOOS 3.0, we also support building inside a Conda environment.  This is the preferred way to build on MacOS.

# Building and Installing LOOS

## For the Impatient


LOOST requires BOOST 1.36 or higher, SCons, and Atlas/LAPACK or other BLAS.
Please refer to the OS-specific instructions below for more details.  For
general advice about configuring LOOS and building in unusual environments, see
the "General Notes" section at the end of this file.

If you are building on a system where the default python is 2.7, you will need to copy custom.py-proto to custom.py, and uncomment the line setting PYTHON_INC (verifying that it's the correct location for your system).

LOOS can then be built using the following command:

    scons

Or installed (to /opt as a default):
    sudo scons install

To install in a user-specified location:

    scons PREFIX=/path/to/install install

To use LOOS, your environment must be first setup:
    (bash)   source /path/to/loos/setup.sh
    (tcsh)   source /path/to/loos/setup.csh



### Build targets

Target | Description
------ | -----------
core   | LOOS Library and PyLOOS
tools  | LOOS Library, Tools, and PyLOOS
all    | LOOS Library, Tools, PyLOOS, and documentation (if necessary), and all Packages (default)
install| Install library, tools, PyLOOS, documentation, and all Packages

### Available Packages (also build targets)

Name    | Description
------  | -----------
ENM     | Elastic Network Models
HBonds  | Hydrogen Bonds Analysis
Conv    | Convergence Analysis
Density | Density/3D Histogram Tools
User    | User-created tools
Python  | PyLOOS scripts


### PyLOOS

The Python interface to LOOS will be included in the build if you have
a recent SWIG (version 2.0 or better) in your standard path as well as
NumPy installed.  Not all operating systems and versions are
supported.  If you need to disable the automatic building of PyLOOS,
use the pyloos flag to scons:

    scons pyloos=0

To build only the core LOOS libraries and PyLOOS, use the following
command:

    scons core


Note that the Optimal Membrane Generator requires PyLOOS.  If you
do not have SWIG installed or disable PyLOOS support, then the OMG
will not be installed.

### Amber NetCDF

LOOS now supports a subset of the Amber NetCDF convention 1.0-B.  Only
coordinates are retrieved and are converted into the default LOOS data
type (i.e. doubles).  Periodic boxes are assumed to be orthogonal and
the angles are currently ignored.

If the netcdf libraries are installed, these will be automatically
detected by SCons and included in the build.  When opening an amber
trajectory file, LOOS will determine if it is a NetCDF file or an
ASCII MDCRD file and act appropriately.

If the netcdf libraries and headers are installed in a non-standard
location, set the NETCDF variable in your custom.py file to point to
the installation.  The specific include and library directories can be
set using the NETCDF_INCLUDE and NETCDF_LIBPATH variables
respectively, and the libraries linked against can be specified using
the NETCDF_LIBS variable.


### Boost

LOOS requires Boost version 1.36 or more recent.  To explicitly
specify an install location, set the BOOST variable in your custom.py
file or on the command line:

    scons BOOST=/usr/local/boost_1_54_0

In some cases, you may need to override either the include directory
or the library directory.  The BOOST_INCLUDE and BOOST_LIBPATH variables
will specify the corresponding directories for the LOOS build.  You
may also explicitly specify which libraries to link against with the
BOOST_LIBS variable.  See custom.py-proto for examples.



### Parallel (multithreaded) ATLAS

If you have a full install of ATLAS including threaded versions of the
BLAS and LAPACK, you can link against these to take advantage of
multiple cores in LOOS.  Copy the custom.py-proto to custom.py and
uncomment/change the appropriate lines.

The MacOS version by default uses the vecLib framework which is
already multithreaded.


### Documentation

However, if you clone LOOS from GitHub, you will need to either:

1. consult the online documentation at http://grossfieldlab.github.io/loos/

2. build a new copy of the documentation.  To do so, you will need to
   install doxygen and graphviz (available in most package managers).

   doxygen


======

# OS Specific Notes

## Conda

Assuming you already have a working install of Anaconda or miniconda, you'll
need to say

    conda create loos
    conda activate loos
    conda install swig scons numpy scipy boost blas libnetcdf

Then, run the newly installed scons to build LOOS

    $CONDA_PREFIX/bin/scons

### Documentation

To build the documentation, you will also require doxygen and graphviz,

    conda install doxygen graphviz
    doxygen


## Fedora

LOOS has been tested on Fedora (64-bit).  We assume you already have
the basic compiler tools installed (i.e. g++).  You will need to
install scons, boost, and atlas:

    sudo dnf install gcc-c++ scons boost-devel atlas-devel netcdf-devel python3-devel swig python3-numpy python3-scipy

LOOS/PyLOOS only supports Python 3, so you must specify which version of numpy
and scipy to use for Fedora 24 and later.  For earlier Fedoras, which don't have
python3 packages for numpy and scipy, you can either install them manually or
use conda.

### Documentation

To build the documentation, you will also require doxygen and graphviz,

    sudo dnf install doxygen graphviz
    doxygen

---

## CentOS

LOOS has been tested with CentOS 6 and 7 (64 bit).  You will need to use the EPEL repo
first,

    sudo yum install epel-release

Then install the packages,

    sudo yum install gcc-c++ scons boost-devel atlas-devel netcdf-devel python-devel swig numpy scipy

### Documentation

To build the documentation, also install:
   sudo yum install doxygen graphviz
   doxygen


### PyLOOS

CentOS 6: Yum install python-devel and pcre-devel, then download
  and install the latest swig.  The version of swig available through
  the package manager is too old.


---

## Ubuntu, Debian, Mint

    sudo apt-get install g++ scons libboost-all-dev libatlas-base-dev libnetcdf-dev swig python3-dev python3-numpy python3-scipy

### Documentation

To build the documentation:
   sudo apt-get install doxygen graphviz
   doxygen

---

## OpenSUSE

As of OpenSUSE 13, there is a pre-build ATLAS package available.  However,
it does not include all of the LAPACK functions LOOS requires.  At this time,
we recommend only installing lapack and blas.  If you install ATLAS, it will
be ignored by the LOOS build.

Using zypper (or your favorite package manager), install the following:

    sudo zypper install gcc-c++ scons boost-devel lapack-devel blas-devel swig netcdf-devel python-numpy python-numpy-devel python-scipy

You should get the blas as a dependency for lapack.  You may also have lapack3
installed by default, however we've found that lapack must also be installed
in order to build LOOS.

### Documentation

To build the documentation:
    sudo zypper install doxygen graphviz
    doxygen

### OpenSUSE 12

The package-manager installed scons is too old.  Download and install SCons
2.0 or better.

---

## MacOS

### IMPORTANT NOTE FOR MACOS 10.11+ USERS !!

There is a problem with using PyLOOS with the new System Integrity
Protection (SIP) enabled (see https://support.apple.com/en-us/HT204899
for more information about SIP).  We are aware of this and working on
a decent solution.  Until then, there are four options for building and
using PyLOOS under El Capitan.

#### Use conda (recommended)

See the conda install instructions above.

#### Build your own Python

The first way is to download and install your own local Python and use
this to run SCons and PyLOOS.  Here again, you have two options: build
and install SCons using that Python, or use the local Python to invoke
the existing SCons.  The latter can be done with the following command
line:
       /path/to/my/own/python `which scons`
Simply putting the new Python in your PATH will not work because SCons
"sanitizes" the PATH before running.  You will also need to be sure to
invoke your PyLOOS scripts with your new local Python.


#### Disable SIP

Although we have not tested this method, it is a common solution
recommended to dealing with unsafe library errors.
Two articles about managing SIP status are:

https://developer.apple.com/library/mac/documentation/Security/Conceptual/System_Integrity_Protection_Guide/ConfiguringSystemIntegrityProtection/ConfiguringSystemIntegrityProtection.html

http://www.macworld.com/article/2986118/security/how-to-modify-system-integrity-protection-in-el-capitan.html


#### Use virtualenv

This is in essence a variant of "install your own python".  If you say
    virtualenv loos-python
it will create a local copy of your python-of-choice in your directory space,
so you can install packages as needed.  Activate that virtual environment using
the activate.csh or activate.sh script in the distribution directory, then use that
to install SCons and LOOS.   You'll then need to make sure you're in that virtualenv
whenever you want to run PyLOOS scripts.

### IMPORTANT NOTE FOR MACOS 10.9 "MAVERICKS" USERS

MacOS 10.9 requires LOOS 2.1 or more recent.  We have also discovered
an incompatibility with Boost installed via Fink in MacOS 10.9.  You
will need to manually download and build a recent version of Boost and
*NOT* use the Fink version.

There is an issue with using Swig and a recent version of MacOS 10.9 that
affects how STL containers are wrapped.  We have disabled the wrapping of
iterator methods for MacOS 10.9 only in order to build PyLOOS.  This means
functions such as begin() and erase() will be unavailable in PyLOOS for
the vectors used in LOOS.

### General Instructions

First, make sure you have the Developer's Tools (i.e. XCode)
installed.  XCode is available for free through the Mac App store.
Next, you will need to install SCons (http://scons.org) and Boost
(http://boost.org) by visiting their websites, downloading the
software, and following their installation instructions.
Alternatively, use fink to install these packages.

#### NetCDF

Download and install the latest hdf5 and netcdf libraries.  If
necessary, set the NETCDF_INCLUDE and NETCDF_LIBPATH variables in your
custom.py file to point to where netcdf is installed.

#### PyLOOS

You will need to download and install a recent version of SWIG first.
If you have installed Boost in a non-standard location, you will need
to make sure that the boost libraries are in your DYLD_LIBRARY_PATH
environment variable.

The default build will use the system Python and Numpy.
Any non-standard locations for python modules can be specified using
the PYTHON_PATH option to scons:
    scons PYTHON_PATH=$HOME/local/lib/python2.7

#### SciPy

Several packages (notably Voronoi) and a few tools (e.g.
cluster-structures.py) depend on scipy.  For recent version of MacOS,
both Scipy should already be installed.  If not, you can always
download it from www.scipy.org.

### Typical Problems


We have seen several instances where LOOS would not build due to
multiple versions of BOOST being installed.  The configuration part of
the build seems to mix components from the different versions
installed.  If your build exits due to errors, verify that you are in
fact using only the BOOST install and libraries you intend (or removed
the excess versions)


---

## Windows/Cygwin (Unsupported)

You will need to install libboost-devel and liblapack-devel along with
the g++ compiler using the cygwin setup.exe program.  (Note: if loos
still does not build, try installing the additional lapack and boost
packages).  There is no atlas package that we could find for cygwin,
so the native lapack/blas will be used instead.

There is no scons install option from setup.exe, so you need to
manually download the latest scons distribution from
http://www.scons.org/download.php (the gzip tar file) and install it.

There is an issue with SCons under cygwin where you must first build
LOOS and -then- install it, i.e.
    scons
    scons PREFIX=/path/to/install/loos install

[Optional]
To include support for Amber/NetCDF trajectory files, install the
libnetcdf-devel library using the cygwin setup.exe installer.


### PyLOOS

Cygwin is NOT SUPPORTED

## Manjaro (Unsupported)

LOOS has been tested with Manjaro 0.8.10.  
Make sure scons, boost, lapack, python, and swig are installed.  Also
install NetCDF, if you want NetCDF support.  LOOS and PyLOOS should build.

---

## Slackware (Unsupported)

LOOS has been tested with Slackware 14.1.  You will need to install,
by whatever means you prefer, lapack, blas, and scons.  LOOS and PyLOOS
should then build.


---

# General Notes


## Customizing the Build

You can override the paths SCons will use for both libraries and
include files by setting the appropriate variables in a "custom.py"
file.  For example, to control where the Boost include files are
located, set the BOOST_INCLUDE variable.

You can also control what libraries are linked against by setting the
appropriate _LIBS variable in your custom.py file.  For example, if
your Boost libraries have a naming convention that the LOOS SConstruct
cannot figure out, you can explicitly set the libraries using the
BOOST_LIBS variable.  These variables take a space-separated list of
library names.  It is important to have *all* required libraries
included in this list.  So for Boost, this would include the regex,
program_options, thread, and system libraries.


If you're using a compiler in a non-standard location (e.g. you have
your own build of the latest and greatest gcc), SCons may not be using
it even though your $PATH is set correctly.  You can force which
compiler is used to build LOOS by setting the CXX variable in your
"custom.py" file.

Note: Settings in the custom.py file can be overridden using the
command-line and the shell environment.  


Sometimes compilers will require multiple environment variables to
work correctly.  In order to handle these cases, the LOOS build will
import all environment variables into SCons before building.  This is
*not* the SCons way to do things, however it makes handling these edge
cases much easier.  In the event that these extra environment
variables cause problems, you can revert to a mostly "clean" build
environment by editing the SConstruct file.  Starting at line 72 with
"env = Environment(...)", uncomment the first invocation and comment
out the second one.


SCons supports building LOOS in parallel.  If you have 4 cores, for
example, use "scons -j4" to use all 4 cores.

Prebuilt documentation for LOOS is provided as part of the
distribution.  This is simply copied into the install directory as
part of installation.  Should you want to build a new version of the
documentation, Doxygen is required.  Moreover, due to an issue we ran
into with SCons, documentation building and installation are
decoupled.  What this means is that you must explicitly build the docs
(i.e. "scons docs") and -then- install, "scons install".  Running
"scons install" will -not- rebuild the documentation, even if it out
of date (or nonexistent).
