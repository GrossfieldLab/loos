# Overview

LOOS is distributed as C++ and python source code, and will need to be compiled
for your particular system.  Both the C++ and python portions have some
external dependencies.

With this in mind, there are 3 routes to installing LOOS:
* conda
* native system Libraries
* roll your own

In each case, you'll still need to compile LOOS on your machine. Conda is the
most portable approach – as far as I can tell, our Conda install works on all
linux and Mac systems; at the moment, *Conda is the only supported approach for
installing LOOS on Mac.*

If you prefer, on Linux you can also use system libraries from your package
manager to satisfy the dependencies.  This file contains the command lines
necessary to install the required libraries for recent versions of Fedora,
Centos, Ubuntu, and OpenSuse.

Finally, if you prefer you can download and build the key dependencies for
yourself, although we can't really help much in that case.

To use LOOS, your environment first be set up.  If you're installing into a
conda environment, you don't have to do anything. Otherwise, you'll need to set
up your environment by sourcing the appropriate setup file. For bash or other
Bourne-like shells, run
```
  source /path/to/loos/setup.sh
```
For csh and tcsh, run
```
source /path/to/loos/setup.csh
```

In addition, you have 2 choices for how to install LOOS. You can either work
with a source tree (use binaries compiled in the directory structure you
downloaded) or an install. Unless you're developing new C++ tools, we suggest an
install is the best choice, particularly if you've built using conda. In that
case, you don't have to do anything to your environment other than activate the
relevant conda environment to use the pre-compiled binaries or write your own
python scripts. If you are developing new tools using C++, or are modifying the LOOS library itself, you'll want to work out of a source tree.

Successful compilation will generate setup.sh and setup.csh in the source
directory, and installation will generate another pair of setup files in the
install directory.  You can have both side by side, and switch from one to the
other by sourcing the appropriate setup file (the scripts prepend the new
directories to the relevant path variables).  _However, if you want to switch
from using a local install to using a conda install, you'll need to say_
```
unset PYTHONPATH
```

_in order to make it work correctly._ If you need to retain other directories in
your PYTHONPATH, you'll have to edit it to remove the directories added by the
setup scripts.

# Compiling using conda for mac or linux

You will need to have a working install of Anaconda or miniconda, available from
https://www.anaconda.com/distribution/

Then, you can run the supplied script to set up a conda environment and build
LOOS

```
   ./conda_build.sh -e loos -j 8 -i
```

This will install packages into an environment loos, creating it if it doesn't
already exist, and will run `scons -j8` (you can supply a different number of
processes if you prefer, eg 2 if you've got a slow machine); the `-i` flag tells
it to do an install. We use conda-forge rather than the default channel, so it's
probably not a great idea to install into an existing environment that uses
other channels. The script will set channel_priority to strict in your
~/.condarc, but you can undo this by removing the following line:

```
channel_priority: strict
```


To create an installation of LOOS, you can say

```
scons install PREFIX=/path/to/loos
```

This defaults to putting LOOS in /opt, but you can choose a different location
either by setting the PREFIX variable, either on the command line or in
custom.py (copy from custom.py-proto to get the idea of what other options are
available).  However, this is not necessary -- you can just as easily work out
of the LOOS source tree, by sourcing `setup.sh` or `setup.csh`, found in the top
level directory.  Installing into the conda distribution is accomplished by
saying

```
scons install PREFIX=$CONDA_PREFIX
```

If you're not using a conda install, we suggest sourcing the appropriate setup
script in your `.bashrc` or `.tcshrc`, so that LOOS will always work for you.

To build the documentation, you will also require doxygen and graphviz,

```
    conda install -c conda-forge doxygen graphviz
    doxygen
```


Going forward, we plan to focus on conda as our preferred environment, and
eventually plan to support direct installation via conda.

Note: if you're updating the source tree of a previous LOOS install, be sure to
remove or rename your `custom.py` file; you shouldn't need it anymore with
conda, and it could mess up the build's search for the correct python, etc.

## Where's my stuff?

Installing into the conda distribution behaves somewhat differently from a
standard install. All executables, both python and c++, get installed into
$CONDA_PREFIX/bin. Python modules, including loos itself, as well as the Voronoi
and OptimalMembraneGenerator packages, are installed into the site-packages
directory of the conda python, eg `$CONDA_PREFIX/lib/python3.8/site-packages/`.
Voronoi and OptimalMembraneGenerator are submodules of loos. Documentation and
examples for OptimalMembraneGenerator are installed in directories inside the
OptimalMembraneGenerator directory, but are probably more easily accessed from
the source distribution.

Building inside a conda environment and installing outside the environment is not a supported configuration.

# Installing using system libraries on supported Linux distributions

The following is a non-exhaustive list of distributions that support LOOS (we've
explicitly tested these systems). It may be possible to adapt the instructions
below to build on older distributions (or you could use conda).

Operating System   | LOOS Support | Notes
----------------   | ------------ | -----
Fedora 29          | yes          |
Fedora 30          | yes          |
Fedora 31          | yes          |
Fedora 32          | yes          |
Ubuntu 16.04 LTS   | yes          | conda-only
Ubuntu 18.04 LTS   | yes          |
Debian 9.9         | yes          |
Centos 7           | yes          | conda-only
Centos 8           | yes          | extra repo
OpenSUSE 15        | yes          |
MacOS X Mojave     | yes          | conda-only
MacOS X Catalina   | yes          | conda-only

You'll need sudo access, or someone who has root access, it order to install
native system libraries. If this isn't possible, your best bet is to install
via conda (see above).

Going forward, we plan to focus on conda as the preferred build environment, as
opposed to using native libraries from the distribution's package manager, and
over time those build instructions may be subject to bit-rot.

After following the instructions specific to your OS, you will need to actually
build LOOS by saying

```
scons -j8
```

"-j8" will run 8 g++ jobs at once, greatly speeding up the build if you've got
a reasonably fast machine.  You can leave this option out, or change the
number, if you prefer.

To create an installation of LOOS, you can say

```
scons install
```

This will install LOOS in /opt. You can also say `scons PREFIX=/path/to/loos
install` to install LOOS in the location of your choosing.

The install will create in $PREFIX the directories bin/, lib/, and include/. It
will also create directories for Voronoi and OptimalMembraneGenerator, but these
are an artifact of the old structuring -- you should use the ones found inside
the lib/ directory instead.

Inside $PREFIX, you will find setup.sh and setup.csh, which will configure your
paths to find this loos distribution -- just source them in your shell.

## Fedora

LOOS has been extensively tested on Fedora.  You will need to install a number
of packages, for instance by using the following command

```
    sudo dnf install gcc-c++ scons boost-devel atlas-devel netcdf-devel python3-devel swig python3-numpy python3-scipy eigen3-devel
```

For versions of Fedora where the default python is python2.7 (anything before
Fedora 31), you will need to copy custom.py-proto to custom.py, and uncomment
the line setting PYTHON_INC (verifying that it's the correct location for your
system), *before* you build LOOS.  Alternatively, you can force scons to run
using python3 by calling it as `/usr/bin/env python3 $(which scons)`.  As of
Fedora 31, this is no longer necessary (the default system python is 3.x).

Fedora versions before 24 don't have python3 packages for numpy and scipy;  to install on these systems, you'll need to either build them manually, or (more efficient) just use conda.

After installing the packages, you can build LOOS using `scons` as described
above.

### Documentation

To build the documentation, you will also require doxygen and graphviz,

```
    sudo dnf install doxygen graphviz
    doxygen
```
---

## CentOS 7

Some combination of the eigen3-devel package and g++ is too old, and leads to
compile errors. To run on CentOS 7, you'll need to build using conda.

## Centos 8

Centos 8 itself doesn't have all of the packages you need to build LOOS -- you need to enable PowerTools and EPEL:

```
sudo yum install https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm
sudo yum config-manager --set-enabled PowerTools
sudo yum install gcc-c++ python3-scons boost-devel atlas-devel netcdf-devel python36 python3-devel swig python3-numpy python3-scipy eigen3-devel

```

For reasons known only to the CentOS team, `scons` is packaged as `scons-3`, so the build and install commands become

```
scons-3
scons-3 install
```


### Documentation

To build the documentation, also install:
```
   sudo yum install doxygen graphviz
   doxygen
```
---

## Ubuntu, Debian, Mint
```
    sudo apt-get install g++ scons libboost-all-dev libboost-regex-dev libatlas-base-dev libnetcdf-dev swig python3-dev python3-numpy python3-scipy libeigen3-dev
```

Copy custom.py-proto to custom.py, and uncomment the line setting PYTHON_INC
(verifying that it's the correct location for your system). Alternatively, you can force scons to run using python3 by saying `/usr/bin/env python3 $(which scons)`.

### Documentation

To build the documentation:
```
     sudo apt-get install doxygen graphviz
     doxygen
```
---

## OpenSUSE

We have tested the build on OpenSuse 15.x
Using zypper (or your favorite package manager), install the following:

```
    sudo zypper install gcc-c++ scons boost-devel lapack-devel blas-devel swig netcdf-devel python-numpy python3-numpy-devel python3-scipy libboost_filesystem1_66_0-devel libboost_program_options1_66_0 libboost_program_options1_66_0-devel libboost_regex1_66_0 libboost_regex1_66_0-devel libboost_system1_66_0-devel libboost_thread1_66_0-devel eigen3-devel
```

You should get the blas as a dependency for lapack.  You may also have lapack3
installed by default, however we've found that lapack must also be installed in
order to build LOOS.

For older distributions, the command line should be very similar, but the
version numbers for boost may differ.

You will need to copy custom.py-proto to custom.py, and uncomment the line
setting PYTHON_INC (verifying that it's the correct location for your system)
*before* running scons. Alternatively, you can force scons to use python3 by
saying `/usr/bin/env python3 $(which scons)`.

### Documentation

To build the documentation:
```
    sudo zypper install doxygen graphviz
    doxygen
```
## General Instructions

This is in case you don't want to use your OS' package manager, and don't want
to use conda.  I'm not sure why anyone would do this, and so these instructions
are really for historical purposes only.

First, make sure you have the Developer's Tools (i.e. XCode, g++, clang, etc)
installed.  Next, you will need to install SCons (http://scons.org) and Boost
(http://boost.org) by visiting their websites, downloading the software, and
following their installation instructions.  

#### NetCDF

Download and install the latest hdf5 and netcdf libraries.  If
necessary, set the NETCDF_INCLUDE and NETCDF_LIBPATH variables in your
custom.py file to point to where netcdf is installed.

#### PyLOOS

You will need to download and install a recent version of SWIG and python3.x
first.  If you have installed Boost in a non-standard location, you will need to
make sure that the boost libraries are in your LD_LIBRARY_PATH (on linux) or
DYLD_LIBRARY_PATH (on OSX) environment variable.

The default build will use the system Python and Numpy.
Any non-standard locations for python modules can be specified using
the PYTHON_INC option to scons:
```
    scons PYTHON_INC=$HOME/local/lib/python3.6
```

You can set this variable by copying  `custom.py-proto` to `custom.py`, and
uncommenting modifying the variables set near the bottom of the file.

#### SciPy

Several packages (notably Voronoi) and a few tools (e.g.
cluster-structures.py) depend on scipy.  You can always
download it from www.scipy.org (and get Numpy while you're at it).

### Typical Problems

We have seen several instances where LOOS would not build due to multiple
versions of BOOST being installed, because the configuration part of the build
seems to mix components from the different versions installed.  If your build
exits due to errors, verify that you are in fact using only the BOOST install
and libraries you intend by examining config.log (or consider removing the
excess versions)

As of LOOS 3.1, this problem should be much rarer, particularly in a Conda
environment, and we'd appreciate hearing about any problems you have.


## Windows (Unsupported)

For Windows 10, your best bet is to use one of the linux subsystems that are installable from Microsoft (e.g. Ubuntu or Debian), then follow the instructions for that linux distribution.  We have anecodotal evidence that this works, but it isn't a supported environment.

It also may be possible to install on windows via conda, but we have not tested this.


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
example, use `scons -j4` to use all 4 cores.

We no longer supply pre-built documentation.  You can either use the docs on the
[GitHub page](http://grossfieldlab.github.io/loos/), or you can build them
yourself.  You'll need to have `doxygen` and `graphviz` installed.  From the
main LOOS directory, run
```
doxygen
```

which will create a new directory `Docs`.  If you open `Docs/html/index.html`,
you'll see an updated version of the docs from the GitHub page (including any
new functions or methods you might have written). The main reason you might want
to do this is if you're adding new classes or methods to the core LOOS library
and want to verify their docs render correctly.


### Build targets

You don't have to build the entirety of LOOS if you don't want to; the
individual components you can choose to build are "build targets", listed
below.  The main reason to do speed up the compile cycle, e.g. if you're
working on a piece of the library, you don't need to rebuild all of the tools
and packages on each compile.

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

Note that the Optimal Membrane Generator and Voronoi packages require PyLOOS.
If you do not have SWIG installed or disable PyLOOS support, then these
packages will not be installed.

### Amber NetCDF

LOOS supports a subset of the Amber NetCDF convention 1.0-B.  Only
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

LOOS definitely requires a Boost distribution more recent than 1.36; we don't
test by manually installing every release, so we can't be sure.  To explicitly
specify an install location if you built it yourself, set the BOOST variable in
your custom.py file or on the command line:

```
    scons BOOST=/usr/local/boost_1_54_0
```

This version of the variable assumes you built Boost yourself, and have it
installed in a single tree (with lib/ and include/ as subdirectories).  If the
tool you used to build Boost broke it up, you may need to override either the
include directory or the library directory.  The BOOST_INCLUDE and BOOST_LIBPATH
variables will specify the corresponding directories for the LOOS build.  You
may also explicitly specify which libraries to link against with the BOOST_LIBS
variable.  See custom.py-proto for examples.

None of this should be necessary if you use Conda.

### Multithreaded linear algebra

If you have a full install of a linear algebra library (e.g. ATLAS, openblas,
etc), you can link against these to take advantage of multiple cores in LOOS.
Copy the `custom.py-proto` to `custom.py` and uncomment/change the appropriate
lines.  Under conda, we use openblas.

### Documentation

You have 2 options for accessing LOOS documentation.  

1. consult the online documentation at http://grossfieldlab.github.io/loos/  
   This is fine if you're not developing new methods for the core library, and if you don't mind needing network access.

2. build a new copy of the documentation.  To do so, you will need to
   install doxygen and graphviz (available in most package managers, including conda).  Then, run

```
   doxygen
```

   from the top-level LOOS directory, and look for the results by accessing
   `Docs/html/index.html`
