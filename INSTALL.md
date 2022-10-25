# Overview


<!--- TODO: This block assumes the conda-forge recipe is finished and working. --->

The best way to install LOOS will depend on how you plan to use it. If you plan to use
current tools and/or write tools in python, by far the easiest way to get it is to install
the package from conda-forge, using the command line

```
    conda install -c conda-forge loos
```

This should install the needed dependencies (netcdf, boost, scipy, sklearn, etc), put the LOOS
binaries and python programs in the conda env's bin, and put the python libraries in the relevant
python's site-packages directory. To use LOOS tools or get the libraries into your environment, you'll need to activate the relevant conda environment:

```
   conda activate ENVNAME
```

where `ENVNAME` is the name of the environment you installed LOOS into. These instructions assume
you've already installed and set up conda; if not you'll need to download either conda or
miniconda from https://docs.conda.io/en/latest/miniconda.html

For the vast majority of users, we expect this is the best approach to take. However, if you plan to
modify the library itself (python or c++) or develop a tool in c++, you'll need a working LOOS tree
to do it. The remainder of the document describes your options for doing so.

There are 3 routes to building LOOS:
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

To use LOOS, your environment needs to be set up.  If you're installing
LOOS into a conda environment, you don't have to do anything.  Otherwise,
you'll need to ensure that the LOOS install location is in your shell's
path (if you are using a non-standard location).


In addition to the conda build, you have 2 choices for how to install LOOS. You can either work
with a build tree or an install. Unless you're developing new C++ tools, we suggest an
install is the best choice, particularly if you've built using conda. In that
case, you don't have to do anything to your environment other than activate the
relevant conda environment to use the pre-compiled binaries or write your own
python scripts. If you are developing new tools using C++, or are modifying the LOOS library itself,
you'll want to work with the build tree.  In this case, you will not only need to set your
shell paths appropriately, but you may also need to set your PYTHONPATH environment
variable to point to the pyloos components.

# Compiling using conda for mac or linux

You will need to have a working install of Anaconda or miniconda, available from
https://www.anaconda.com/distribution/

Then, you can run the supplied script to set up a conda environment and build
LOOS

```
   ./conda_build.sh -e loos -j 8 -i
```

This will install packages into an environment loos, creating it if it doesn't
already exist, and will run cmake on 8 cores (you can supply a different number of
processes if you prefer, eg 2 if you've got a slow machine); the `-i` flag tells
it to do an install. We use conda-forge rather than the default channel, so it's
probably not a great idea to install into an existing environment that uses
other channels. The script will set channel_priority to strict in your
~/.condarc, but you can undo this by removing the following line:

```
channel_priority: strict
```


If you want to build everything by
hand, but with conda, first set up your conda environment,

```
conda create -n loos -c conda-forge python=3 swig=4 cmake numpy scipy scikit-learn boost openblas libnetcdf=4.8.1 lapack compilers eigen
conda activate loos
```

Then build loos with the following.  You can use this same procedure
if you make a change to your loos distribution (or run `git pull`) and want to
rebuild it (note: you can add `-j n` to build with `n` processes)

```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX ..
cmake --build .
cmake --install .
```

If you are using a conda environment, be sure to activate it _before_ running
any `cmake` commands.  If you don't specify an install location, the default
will be `/usr/local/bin`. To install outside the conda environment, you should
replace `$CONDA_PREFIX` in the command above  with your preferred target. cmake
will install the python library  components in the python tree of the python
used in the environment running cmake.

To build the documentation, you will also require doxygen and graphviz,

```
conda install -c conda-forge doxygen graphviz

cmake ..
cmake --build . --target=docs
```

Please note that the documentation is currently only supported in the build and is not installed.

If you want to configure how loos is built, you can either set the `cmake` variables directly,
or use the graphical configuration tool `ccmake`, e.g.

```
cd build
ccmake ..
```


Going forward, we plan to focus on conda as our preferred environment, and
eventually plan to support direct installation via conda.

If you want to clean out the build generated by cmake, cd into the build directory and 
say 
```
cmake --build . --target clean
```

## Upgrading

If you updated an existing git repo from loos 3.x, we suggest running `scons -c;
scons -c config` before pulling in changes from greater than or equal to 4.0.0. We also suggest building into a fresh conda environment. 

However, we know that sometimes this is a pain, and some have had success in
retaining their envs.  If you are trying to rebuild the new loos in an extant
conda env, you should still run the aforementioned scons commands to remove
kruft from your 3x install before adding anything from 4.x. This applies even if
you're switching from a local build to the new conda-forge binaries if you are
trying to re-use your env.

If you accidentally installed the cmake build overtop of an older loos build in the same env, one of us managed to get to a working environment by:
1. Following the 4.x uninstall instructions above.
2. Checking out tagged release 3.3, then running `scons -c; scons -c config`.
3. Manually inspecting the contents of the `/path/to/miniconda3/envs/my-old-loos-env` to make sure this 'got' everything (for example, using `find` to look for files with `loos` or `scons` somewhere in the name). There may not be any, but an ounce of prevention...
4. Returning to the loos repo, checking out the version you wanted to build in the first place, and proceeding with the cmake install.

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

## Uninstalling

If you want to remove loos from your environment, you can run,

```
pip uninstall loos
```

## Upgrading

If you are upgrading from a 3.x version, and you would like to retain your conda env, follow the uninstall instructions for your current `scons` built version:

## Working with the build (i.e. not installing)

If you want to work with the build directory as your loos installation, rather than actually installing it,
you can find all of the binary executables either under `Tools` or in a subdirectory of `Packages`.  Note
that the PyLOOS Packages do not have anything that has to be "built", so those stay in the source tree
until it's time to install.

To use pyloos, you will want to set your `$PYTHONPATH` environment variable to
point into the appropriate build directory,

```
export PYTHONPATH=/path/to/build/src/pyloos/src
```

## Alternative build systems

Using build systems other than make is possible with the new CMake version of loos.  For example,
if you'd like to use Ninja instead, you can install it using conda,

```
conda activate loos
conda install ninja
```

And then tell cmake to use it (but be sure to clear out your build directory first),

```
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX .. -GNinja
cmake --build . -j4
cmake --install .
```


# Installing using system libraries on supported Linux distributions

The following is a non-exhaustive list of distributions that support LOOS (we've
explicitly tested these systems). It may be possible to adapt the instructions
below to build on older distributions (or you could use conda).

Operating System   | LOOS Support | Notes
----------------   | ------------ | -----
Fedora 36          | yes          | 
Ubuntu 20.04 LTS   | yes          | conda-only
Debian 10.x        | yes          | conda-only
Debian 11          | yes          | 
Centos 8           | yes          | extra repo
OpenSUSE 15        | yes          | *
MacOS X Mojave     | yes          | conda-only
MacOS X Catalina   | yes          | conda-only
MacOS X Monterrey  | yes          | conda-only

You'll need sudo access, or someone who has root access, it order to install
native system libraries. If this isn't possible, your best bet is to install
via conda (see above).

Going forward, we plan to focus on conda as the preferred build environment, as
opposed to using native libraries from the distribution's package manager, and
over time those build instructions may be subject to bit-rot.

After following the instructions specific to your OS, you will need to actually
build LOOS by saying

```
cd /path/to/loos/source/distribution
mkdir build
cd build
cmake ..
cmake --build . -j8
```

"-j8" will run 8 g++ jobs at once, greatly speeding up the build if you've got
a reasonably fast machine.  You can leave this option out, or change the
number, if you prefer.

To create an installation of LOOS, you can say

```
cmake --install .
```

This will install LOOS in `/usr/local`.  To specify where to install LOOS, use the following,
```
cmake -DCMAKE_INSTALL_PREFIX=/path/to/loos/install ..
cmake --install . -j8
```

Note that all PyLOOS and python components will be installed to your current python, in the site-packages directory.  
We recommend you either create a conda environment and use that python, or a python virtual environment.


## Fedora

LOOS has been tested on Fedora.  You will need to install a number
of packages, for instance by using the following command

```
    sudo dnf install cmake gcc-c++  boost-devel atlas-devel netcdf-devel python3-devel python3-pip swig python3-numpy python3-scipy eigen3-devel python3-scikit-learn
```

After installing the packages, you can build LOOS using `cmake` as described
above.

### Documentation

To build the documentation, you will also require doxygen and graphviz,

```
    sudo dnf install doxygen graphviz
```
---

## CentOS 7

Some combination of the eigen3-devel package and g++ is too old, and leads to
compile errors. To run on CentOS 7, you'll need to build using conda.

## Centos 8

Centos 8.5 itself doesn't have all of the packages you need to build LOOS -- you need to enable PowerTools and EPEL:

```
sudo yum install https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm
yum install dnf-plugins-core
sudo yum config-manager --set-enabled PowerTools
sudo yum install cmake gcc-c++ boost-devel atlas-devel netcdf-devel python36 python3-devel swig python3-numpy python3-scipy eigen3-devel
```

Note: for reasons I don't understand, sometimes the `yum config-manager` line wants it
written as "powertools" with no capitals. If you get a message saying the repo isn't there,
try that.

Note: Centos8 doesn't have a python3-scikit-learn package, so a couple of LOOS tools
that depend on it won't work. If you need scikit-learn, you'll have to build under conda instead.



### Documentation

To build the documentation, also install:
```
   sudo yum install doxygen graphviz
```
---

## Ubuntu, Debian, Mint
```
    sudo apt-get install cmake g++ libboost-all-dev libboost-regex-dev libatlas-base-dev libnetcdf-dev swig python3-dev python3-pip python3-numpy python3-scipy libeigen3-dev python3-sklearn python3-sklearn-lib
```


### Documentation

To build the documentation:
```
     sudo apt-get install doxygen graphviz
     doxygen
```
---

## OpenSUSE

We have tested the build on OpenSuse 15.5
Using zypper (or your favorite package manager), install the following:

```
    sudo zypper install cmake gcc-c++ boost-devel lapack-devel blas-devel swig netcdf-devel python3-pip python3-numpy-devel python3-scipy python3-scikit-learn libboost_filesystem1_66_0-devel libboost_program_options1_66_0 libboost_program_options1_66_0-devel libboost_regex1_66_0 libboost_regex1_66_0-devel libboost_system1_66_0-devel libboost_thread1_66_0-devel eigen3-devel
```

NOTE: at least some versions of OpenSuse15 have broken numpy packages for python 3.
In this case, trying to import loos (or just numpy itself) will give a missing
symbol error. In this case, you'll need to go with a conda build.

You should get the blas as a dependency for lapack.  You may also have lapack3
installed by default, however we've found that lapack must also be installed in
order to build LOOS.

For older distributions, the command line should be very similar, but the
version numbers for boost may differ and cmake might not be new enough. Or,
just use conda.


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
installed.  Next, you will need to install CMake and Boost
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

The default build will attempt to find the correct python and numpy.

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
install, by whatever means you prefer, lapack, blas, and cmake.  LOOS and PyLOOS
should then build.

Or, just use conda.


---

# General Notes


## Customizing the Build



We no longer supply pre-built documentation.  You can either use the docs on the
[GitHub page](http://grossfieldlab.github.io/loos/), or you can build them
yourself.  Follow the instructions above for building the documentation.
You will find it in the `build/html` directory.  If you open `build/html/index.html`,
you'll see an updated version of the docs from the GitHub page (including any
new functions or methods you might have written). The main reason you might want
to do this is if you're adding new classes or methods to the core LOOS library
and want to verify their docs render correctly.

### Amber NetCDF

LOOS supports a subset of the Amber NetCDF convention 1.0-B.  Only
coordinates are retrieved and are converted into the default LOOS data
type (i.e. doubles).  Periodic boxes are assumed to be orthogonal and
the angles are currently ignored.

If the netcdf libraries are installed, these will be automatically
detected by CMake and included in the build.  When opening an amber
trajectory file, LOOS will determine if it is a NetCDF file or an
ASCII MDCRD file and act appropriately.

If the netcdf libraries and headers are installed in a non-standard
location, you can set the `NETCDF_INCLUDES` and `NETCDF_LIBRARIES` CMake variables.


### Boost

LOOS definitely requires a Boost distribution more recent than 1.36; we don't
test by manually installing every release, so we can't be sure.  If you need
to manually specify where Boost has been installed, you can use the following
CMake variables: `BOOST_ROOT` or `BOOST_INCLUDEDIR` and `BOOST_LIBRARYDIR`.

None of this should be necessary if you use Conda.

### Multithreaded linear algebra

If you need to customize the blas/lapack used, you can set the `BLA_VENDOR` variable
(see [FindLAPACK]https://cmake.org/cmake/help/latest/module/FindLAPACK.html)

Under conda, we use openblas.

### Documentation

You have 2 options for accessing LOOS documentation.  

1. consult the online documentation at http://grossfieldlab.github.io/loos/  
   This is fine if you're not developing new methods for the core library, and if you don't mind needing network access.

2. build a new copy of the documentation.  To do so, you will need to
   install doxygen and graphviz (available in most package managers, including conda).  Then, build as described above.

   Look for the results by accessing
`build/html/index.html`
