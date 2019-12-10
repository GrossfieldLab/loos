# Operating System Compatibility

                   |     
Operating System   | LOOS Support | PyLOOS Support | Notes
----------------   | ------------ | -------------- | -----
Fedora 18          | yes          | yes            | Deprecated
Fedora 19          | yes          | yes            | Deprecated
Fedora 20          | yes          | yes            | Deprecated
Fedora 21          | yes          | yes            | Deprecated
Fedora 22          | yes          | yes            | Deprecated
Fedora 23          | yes          | yes            | Deprecated
Fedora 24          | yes          | yes            | Deprecated
Fedora 25          | yes          | yes            | Deprecated
Fedora 26          | yes          | yes            | Deprecated
Fedora 27          | yes          | yes            | Deprecated
Fedora 28          | yes          | yes            |
Fedora 29          | yes          | yes            |
Fedora 30          | yes          | yes            |
Fedora 31          | yes          | yes            |
Ubuntu 12.04 LTS   | yes          | yes            | Deprecated
Ubuntu 14.04 LTS   | yes          | yes            | Deprecated
Ubuntu 15.04       | yes          | yes            | Deprecated
Ubuntu 15.10       | yes          | yes            | Deprecated
Ubuntu 16.04 LTS   | yes          | yes            | Deprecated
Ubuntu 18.04       | yes          | yes            |
Debian 7.8         | yes          | yes            | Deprecated
Debian 8.1         | yes          | yes            | Deprecated
Debian 9.8         | yes          | yes            |
Centos 7           | yes          | yes            |
OpenSUSE 12        | yes          | yes            | Deprecated
OpenSUSE 13        | yes          | yes            | Deprecated
OpenSUSE 15        | yes          | yes            |
MacOS X            | yes          | yes            | See OS notes


* Deprecated: We used to support this configuration, but no longer test it.  It may still work. Your best bet is to use conda.

As of LOOS 3.0, we also support building inside a Conda environment.  This is the preferred way to build on MacOS, and on any Linux environment that is not supported.  LOOS 3.1 features an extensively reworked build system, which uses conda much more natively and should be far more robust.

Going forward, we plan to focus on conda as the preferred build environment, as opposed to using native libraries from the distribution's package manager, and over time those build instructions may be subject to bit-rot.


# Building and Installing LOOS

## For the Impatient

For the really impatient, you can just run `conda_build.sh`, found in the main LOOS directory, which will create an appropriate conda environment and build LOOS in it for you.

LOOS requires BOOST 1.36 or higher, SCons, and Atlas/LAPACK or other BLAS.
Please refer to the OS-specific instructions below for more details.  For
general advice about configuring LOOS and building in unusual environments, see
the "General Notes" section at the end of this file.

If you are building on a system where the default python is 2.7, you will need to copy custom.py-proto to custom.py, and uncomment the line setting PYTHON_INC (verifying that it's the correct location for your system).

LOOS can then be built using the following command:
```
    scons
```

Or installed (to /opt as a default):
```
    sudo scons install
```

To install in a user-specified location:

```
    scons PREFIX=/path/to/install install
```

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
NumPy installed.  If you need to disable the automatic building of PyLOOS,
use the pyloos flag to scons:

```
    scons pyloos=0
```

To build only the core LOOS libraries and PyLOOS, use the following
command:

```
    scons core
```

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

```
    scons BOOST=/usr/local/boost_1_54_0
```

In some cases, you may need to override either the include directory
or the library directory.  The BOOST_INCLUDE and BOOST_LIBPATH variables
will specify the corresponding directories for the LOOS build.  You
may also explicitly specify which libraries to link against with the
BOOST_LIBS variable.  See custom.py-proto for examples.

None of this should be necessary if you use Conda.



### Parallel (multithreaded) ATLAS

If you have a full install of ATLAS including threaded versions of the
BLAS and LAPACK, you can link against these to take advantage of
multiple cores in LOOS.  Copy the custom.py-proto to custom.py and
uncomment/change the appropriate lines.



### Documentation

You have 2 options for accessing LOOS documentation.  

1. consult the online documentation at http://grossfieldlab.github.io/loos/  
   This is fine if you're not developing new methods for the core library, and if you don't mind needing network access.

2. build a new copy of the documentation.  To do so, you will need to
   install doxygen and graphviz (available in most package managers).  Then, run

```
   doxygen
```

   from the top-level LOOS directory, and look for the results by accessing
   `Docs/html/index.html`


======

# OS Specific Notes

## Conda

You will need to have a working install of Anaconda or miniconda, available from
https://www.anaconda.com/distribution/

Then, you can run the supplied script to set up a conda environment and build LOOS

```
   ./conda_build.sh loos 8
```

This will install packages into an environment loos, creating it if it doesn't
already exist, and will run `scons -j8` (you can supply a different number of
processes if you prefer, eg 2 if you've got a slow machine).  We use conda-forge
rather than the default channel, so it's probably not a great idea to install
into an existing environment that uses other channels. The script will set
channel_priority to strict in your ~/.condarc, but you can undo this by removing
the following line:

```
channel_priority: strict
```

As of version 3.1, we've significantly redone the build scripts so this should work robustly.  However, if the build fails to find something (or you want to use a version of a library from outside of conda), you can copy custom.py-proto to custom.py, and set the variables (e.g. BOOST, NETCDF) to point to the locations of the libraries you want to use.

To create an installation of LOOS, you can say

```
scons install
```

This defaults to putting LOOS in /opt, but you can choose a different location either by setting the PREFIX variable, either on the command line or in custom.py.  However, this is not necessary -- you can just as easily work out of the LOOS source tree, by sourcing `setup.sh` or `setup.csh`, found in the top level directory.


### Documentation

To build the documentation, you will also require doxygen and graphviz,

```
    conda install doxygen graphviz
    doxygen
```

## Fedora

LOOS has been tested on Fedora (64-bit).  We assume you already have
the basic compiler tools installed (i.e. g++).  You will need to
install scons, boost, and atlas:

```
    sudo dnf install gcc-c++ scons boost-devel atlas-devel netcdf-devel python3-devel swig python3-numpy python3-scipy
```

You may need to copy custom.py-proto to custom.py, and uncomment the line
setting PYTHON_INC (verifying that it's the correct location for your system).
As of Fedora 31, this is no longer necessary (the default system python is
3.x).

LOOS/PyLOOS only supports Python 3, so you must specify which version of numpy
and scipy to use for Fedora 24 and later.  For earlier Fedoras, which don't have
python3 packages for numpy and scipy, you can either install them manually or
use conda.

### Documentation

To build the documentation, you will also require doxygen and graphviz,

```
    sudo dnf install doxygen graphviz
    doxygen
```
---

## CentOS 7

You'll need the epel repository in order to get python3 versions of numpy and scipy:

```
    yum install https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm
```

Then install the packages

```
    sudo yum install gcc-c++ scons boost-devel atlas-devel netcdf-devel python36 python36-devel swig python36-numpy python36-scipy
```

You may need to copy custom.py-proto to custom.py, and uncomment the line setting PYTHON_INC (verifying that it's the correct location for your system).


### Documentation

To build the documentation, also install:
```
   sudo yum install doxygen graphviz
   doxygen
```
---

## Ubuntu, Debian, Mint
```
    sudo apt-get install g++ scons libboost-all-dev libatlas-base-dev libnetcdf-dev swig python3-dev python3-numpy python3-scipy
```

Copy custom.py-proto to custom.py, and uncomment the line setting PYTHON_INC (verifying that it's the correct location for your system).

### Documentation

To build the documentation:
```
     sudo apt-get install doxygen graphviz
     doxygen
```
---

## OpenSUSE

As of OpenSUSE 13, there is a pre-build ATLAS package available.  However,
it does not include all of the LAPACK functions LOOS requires.  At this time,
we recommend only installing lapack and blas.  If you install ATLAS, it will
be ignored by the LOOS build.

Using zypper (or your favorite package manager), install the following:

```
    sudo zypper install gcc-c++ scons boost-devel lapack-devel blas-devel swig netcdf-devel python-numpy python3-numpy-devel python3-scipy libboost_filesystem1_66_0-devel libboost_program_options1_66_0 libboost_program_options1_66_0-devel libboost_regex1_66_0 libboost_regex1_66_0-dev libboost_system1_66_0-devel libboost_thread1_66_0-devel
```

You should get the blas as a dependency for lapack.  You may also have lapack3
installed by default, however we've found that lapack must also be installed
in order to build LOOS.

You may need to copy custom.py-proto to custom.py, and uncomment the line setting PYTHON_INC (verifying that it's the correct location for your system).

### Documentation

To build the documentation:
```
    sudo zypper install doxygen graphviz
    doxygen
```
### OpenSUSE 12

The package-manager installed scons is too old.  Download and install SCons
2.0 or better.

---

## MacOS

We only support OS X via conda -- see the conda install instructions above.  We have worked to remove external dependencies, so it shouldn't be necessary to have XCode installed.


### General Instructions

This is in case you don't want to use you OS' package manager, and don't want to use conda.  I'm not sure why anyone would do this, and so these instructions are really for historical purposes only.

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

We have seen several instances where LOOS would not build due to multiple
versions of BOOST being installed, because the configuration part of the build
seems to mix components from the different versions installed.  If your build
exits due to errors, verify that you are in fact using only the BOOST install
and libraries you intend by examining config.log (or consider removing the
excess versions)

As of LOOS 3.1, this problem should be much rarer, particularly in a Conda
environment, and we'd appreciate hearing about any problems you have.


---

## Windows (Unsupported)

For Windows 10, your best bet is to use one of the linux subsystems that are installable from Microsoft (e.g. Ubuntu or Debian), then follow the instructions for that linux distribution.  We have anecodotal evidence that this works, but it isn't a supported environment.

It also may be possible to install on windows via conda, but we have not tested this.

---

## Slackware (Unsupported)

Older versions of LOOS have been tested with Slackware 14.1.  You will need to
install, by whatever means you prefer, lapack, blas, and scons.  LOOS and PyLOOS
should then build.

Or, just use conda.


---

# General Notes


## Customizing the Build

You can override the paths SCons will use for both libraries and
include files by setting the appropriate variables in a "custom.py"
file.  For example, to control where the Boost include files are
located, set the BOOST_INCLUDE variable.

You can also control what libraries are linked against by setting the
appropriate `_LIBS` variable in your custom.py file.  For example, if
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

SCons supports building LOOS in parallel.  If you have 4 cores, for
example, use "scons -j4" to use all 4 cores.

We no longer supply pre-built documentation.  You can either use the docs on the
[GitHub page](http://grossfieldlab.github.io/loos/), or you can build them
yourself.  You'll need to have `doxygen` and `graphviz` installed.  From the
main LOOS directory, run
```
doxygen
```

which will create a new directory `Docs`.  If you open `Docs/html/index.html`,
you'll see an updated version of the docs from the GitHub page (including any
new functions or methods you might have written).
