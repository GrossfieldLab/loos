# The Lightweight Object Oriented Structure analysis library

![Build Status](https://github.com/GrossfieldLab/loos/workflows/CI/badge.svg)

## Getting Help
* Raise an [issue](https://github.com/GrossfieldLab/loos/issues) here on GitHub.
* Ask a question on the GitHub [discussion board](https://github.com/GrossfieldLab/loos/discussions)
* Look at the [contributions guide](https://github.com/GrossfieldLab/loos/blob/main/.github/CONTRIBUTING.md)
* You can also contact us directly at loos.maintainer@gmail.com

## Documentation


Documentation is currently available online at
http://grossfieldlab.github.io/loos/.

These pages include brief summaries of most of the tools available with LOOS, as
well as auto-generated class and method level documentation for developers. If
you want a local copy of the documentation, you can build it by running doxygen
in the main LOOS directory.  See the [INSTALL](INSTALL.md) file for more
details.

Additional documentation is available on the [GitHub wiki](https://github.com/GrossfieldLab/loos/wiki), including [slides from a talk introducing LOOS](http://membrane.urmc.rochester.edu/sites/default/files/loos.pdf) and brief articles focused on [how to use LOOS](https://github.com/GrossfieldLab/loos/wiki/Tutorials-for-Users) and how to [develop with LOOS](https://github.com/GrossfieldLab/loos/wiki/Tutorials-for-Developers).

## Referencing loos

If you use LOOS in published work, please cite these two papers:

- Romo, T.D., Grossfield, A. "LOOS: An extensible platform for the structural analysis of simulations." 31st Annual International Conference of the IEEE EMBS (2009): 2332-2335

- Romo, T. D., Leioatts, N. and Grossfield, A., "Lightweight Object Oriented Structure Analysis: Tools for building tools to analyze molecular dynamics simulations", J. Comput. Chem. (2014): 2305-2318


## Using LOOS

Welcome to the GitHub version of LOOS.  This repository was
created by converting our old SVN repository.  The SVN feature
branches exist as refs/tags, and appear in GitHub as releases,
but they should not be confused with the actual LOOS releases.  Those
can be found using semantic versioning (e.g. release-2.3.1,
release-2.3.0, ...)

We don't release all that often, but we maintain the main branch in a correct
and usable state.  All development is done in other branches, and merged once we
believe it's correct.


For help with installing LOOS, please see the [INSTALL.md](INSTALL.md) file.  For
more details about what has changed in LOOS, see the [ChangeLog](ChangeLog) file.

## Layout of this repository

The structure of this repository changed substantively for the 4.0 release. LOOS is 
simultaneously a C++ library for analyzing molecular dynamics simulations, a python 
wrapper for that library, and a suite of tools in C++ and python. The main directory 
contains a src/ containing the code for the C++ library, Tools/ containing applications
written in C++, and Packages/, which contains several sets of tools with extra code 
that don't really fit inside the library. These include Packages/Voronoi, which performs
voronoi analysis on membrane systems, Packages/Clustering, which does clustering of 
macromolecules (although there are other clustering tools found elsewhere in LOOS), 
Packages/Convergence, which focuses on estimating uncertainty and statistical error, 
Packages/ElasticNetworks, which implements a variety of elastic network analysis, 
Packages/HydrogenBonds, which does hydrogen bond analysis, Packages/PyLOOS, which contains 
tools written in Python, and Packages/User, is a location for users to easily 
add new C++ programs. Code for implementing the python bindings (including library
code associated with packages) is found in src/loos. Finally, the top-level share/ directory
contains several data files associated with particular tools (e.g. suite definitions for the
`rna_suites` tool. 

Users looking for simple example programs on which they can base new code can start with 
Packages/User (for C++ tools) and Packages/PyLOOS (for Python tools). If you are considering 
writing a new tool, we suggest starting with the python. Obviously, writing python is easier,
and generally speaking the performance hit for python vs. C++ isn't prohibitive. You can find
resources for developing with LOOS on the [GitHub wiki](https://github.com/GrossfieldLab/loos/wiki),
particularly the [Tutorials for Developers](https://github.com/GrossfieldLab/loos/wiki/Tutorials-for-Developers).

### Release 4.0.3

LOOS 4.0 is now officially "in production" and can be installable either via conda-forge or
by building from source.  The first few versions had minor issues either related to
the new build system (cmake) or the conda-forge packaging.  The current release is 4.0.3
and we recommend that everyone update to this as soon as possible.

If you have been using LOOS prior to 4.0, particulary with the conda build system, then the
safest way to upgrade is to recreate your conda environment and install the new conda-forge
package.  If you chose to build from source, but had previously installed LOOS into your conda
environment, then we still recommend that you delete and recreate your environment, _then_ build
and install the latest LOOS.

If you have installed a previous version of LOOS from conda-forge, then you can update to
the latest with,

```
conda activate loos
conda update loos
```



### Release pre-4.0.0

This release involves a complete switch from using SCons for building LOOS to CMake.
In general, you should find the build process simpler and configuration easier.

This is a pre-release version, so it should be considered experimental and not used
for production.  Features may change without notice.

### Release 3.3.0

This release has a new and improved facility for reweighting (from Louis Smith), methods for calculating logistic and hard-cutoff lateral densities in AtomicGroup, several new tools for working with clusters in the Clustering package, and dipole computation in membrane_map. Also, a new tool (ocf) for calculating a quantity related to persistence length (also from Louis Smith).

There are also a number of bugfixes:
 - Incorrect bonds removed from add_molecule.py
 - Fixed units for periodic box read by gmxdump2pdb.pl
 - Update to handle new constraints format in gromacs 2020.3+ in gmxdump2pdb.pl
 - Fixed charge units when reading Amber prmtop files
 - OMG correctly handles 4- and 5-site water models
 - In subsetter and merge-traj, the --regex and --scanf flags now imply --sort

 There were also a few housekeeping updates. We migrated continuous integration from Travis-CI to GitHub Actions, Louis Smith contributed a *much* cleaner .gitignore file, and we fixed a potential build-bug where the C++ file generated by swig had different names depending on which compiler/swig version was used. We also added scikit-learn to the list of conda packages installed by conda_build.sh, and fixed a silly bug in its detection of the available compilers.

### Release 3.1.0

This is another major release. The biggest change is a totally reworked build
system, which should make the process of installing LOOS much easier.  In
particular, we've greatly improved our support for installing under Conda.  See
[INSTALL.md](INSTALL.md) for detailed instructions.

We've added a new package, Clustering. At the moment, it only has one tool (a fast k-medioids tool), but several others are planned.  There are also a number of improvements to other tools.

### Release 3.0.0

This is a major release.  Most notably, from the development side we've converted from python 2.7 to python 3.x for better long-term interoperability. We've also finally dropped SourceForge as a distribution platform (the old site forwards here), and are deprecating the old mailing lists.  Instead, to keep abreast of LOOS development, we suggest following the project on GitHub.  Similarly, we recommend raising issues on GitHub as the best way to ask for support, although emailing us directly at loos.maintainer@gmail.com will also work.


### Release 2.3.3 (beta)

Note: Checking out a beta release means you are checking out a version
of LOOS that is under active development.  This may include build issues,
tools not working, and undocumented "features."

Major changes in LOOS include better DCD handling, support for multiple
trajectories in some tools (and at the API level), as well as a new parser
for specifying frame ranges for tools.

LOOS now has the ability to handle DCD trajectories with
a 0 frame count in the header (fixdcd is no longer required for this
case).  The count will be estimated based on the model-size and trajectory
file size.

A new subclass of Trajectory has been added called MultiTrajectory.  A
MultiTrajectory object may contain multiple pTraj's and treats them as
one giant trajectory.  Each sub-trajectory can have its own skip and stride.
In addition, there is a MultiTrajOptions class for handling multiple
trajectories in a tool.

The parse for specifying ranges has been upgraded to use Boost Spirit.  This
should make it more robust.  In addition, you now no longer need to know
the length of a trajectory to use a range.  An empty "field" will be
filled in with the appropriate value.  For example, to skip the first 10 frames,
then take every other frame until the end of the trajectory, use:
     toolname -r 10:2: model.psf trajectory.dcd

NOTE: The short options to subsetter have been changed to be consistent with the
new MultiTrajOptions set of options: "-i" is used for stride and "-k" for skip.
The long options have not changed.

Additional changes to LOOS include the addition of a new
lipid_survival tool and a multi-rmsds tool.  Proper (full) support for
atom inequalities in Python has been addressed.  A new reimaging mode
has been added to subsetter (--reimage=zealous) that fixes some
issues coming from Gromacs.  The output of dibmops has been changed to
have a "0" in bins with no data rather than "-1".  Finally, a number of
bugs have been fixed.  See the ChangeLog for more details.


### RELEASE 2.3.2

This release includes a number of changes to support our migration to
GitHub as well as some important bug fixes and additions.

We have reorganized the LOOS source code so that the core library now
resides in the "src" directory.  The shared library that's built (along
with the python code) is still copied to the top-level LOOS directory.

The Doxygen-based documentation is now handled a little bit differently.
When cloning from GitHub, the documentation will be automatically built
by SCons.  This means doxygen and graphviz are now required to build LOOS.
If you download a release from GitHub, then you can also download the pre-built
tar file containing the documentation.  If this is in the top LOOS
source directory, then SCons will see it and unpack it for you.  Finally,
if you download a release from SourceForge, the the pre-built documentation
is bundled with the release.  In all cases, you can always find the docs
for the current release online at http://grossfieldlab.github.io/loos/
or http://loos.sourceforge.net

Bugs fixed in this release include some important fixes to the OMG,
a bug affecting aligner and membrane_map when aligning to a reference
structure, making XTCWriter available to PyLOOS, and a few more minor
fixes.  See the [INSTALL](INSTALL.md) file for more details.

New tools in this release include verap, a quick vertical area profile
tool, cylindrical-thickness, and inside_helices.  New features include
providing support for manually mapping molecule names to segids in the
gmxdump2pdb tool and support for writing GRO files.


### RELEASE 2.3.1

This release of LOOS is a bug-fix release.  A bug was discovered in
membrane_map that caused the z-axis of all coordinates to be set to
0.  This affects height and vector calculations, but not density.

The build system has been changed slightly.  Some compilers require
multiple environment variables to work correctly.  In order to handle
these cases, the LOOS build will import all environment variables into
SCons before building.  This is *not* the SCons way to do things,
however it makes handling these edge cases much easier.  In the event
that these extra environment variables cause problems, you can revert
to a mostly "clean" build environment by editing the SConstruct file.
Starting at line 72 with "env = Environment(...)", uncomment the first
invocation and comment out the second one.

In addition, improvements were made to the documentation.




### RELEASE 2.3.0

This release of LOOS comes with a number of major changes, most
related to PyLOOS.  A detailed description of these changes is
available in the new "what's new in PyLOOS" section of the
documentation, and is listed in the "Quick Links" off the main page.

NUMPY IS NOW REQUIRED FOR BUILDING PYLOOS.  For Linux users, this will
just be a package install.  For Mac users, the system numpy will be
used.  See the [INSTALL](INSTALL.md) file for more information.

PyLOOS has been reorganized and new trajectory handling classes have
been added.  All Python-based components of PyLOOS are now in the
loos.pyloos module.  The new trajectory classes are easier to use than
the old PyTraj classes (they no longer require you to explicitly wrap
a loos.Trajectory object).  In addition, you can combine multiple
trajectories into a large "virtual" trajectory.  The frames of the
composite trajectory can be iteratively aligned (a la aligner) or
aligned to a reference structure.  See the documentation for pyloos
under the Namespace tab or search for Trajectory.

Additional PyLOOS changes include bug fixes for the iterativeAlignment
function.  The splitting functions in AtomicGroup now return a Python
list rather than an AtomicGroupVector.  Coordinates can be extracted
from an AtomicGroup or a Trajectory as a NumPy matrix.  An SVD/PCA
function has also been added to PyLOOS.  All old PyLOOS functions that
begin or end with Py are now deprecated and will be removed in the
next release.  These include PyTraj, PyAlignedTraj, and
iterativeAlignmentPy.  Finally, a k-means clustering program (using
PyLOOS) has been added.  Also, we've added a new package (Voronoi), which
will help users analyze membrane/membrane protein simulations.

There are a number of significant non-PyLOOS changes as well.  An
"index" keyword has been added to pick atoms based on their position
in a model (rather than atom id).  The rmsds tool (all-to-all RMSD)
has been significantly optimized and can now be run multithreaded.
The default is to use only a single thread.  If you built LOOS with a
multithreaded math library, be aware of possible conflicts (though
that shouldn't happen with rmsds).  Last, but not least, a number of
minor bug fixes are included.  See the ChangeLog file or the change
log in the documentation for more details.
