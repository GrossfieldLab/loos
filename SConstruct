#!/usr/bin/env python3
#  This file is part of LOOS.
#
#  LOOS (Lightweight Object-Oriented Structure library)
#  Copyright (c) 2008-2009, Tod D. Romo
#  Department of Biochemistry and Biophysics
#  School of Medicine & Dentistry, University of Rochester
#
#  This package (LOOS) is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation under version 3 of the License.
#
#  This package is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import os
import os.path
import glob
import platform
import re
import subprocess
from time import strftime
import shutil
import distutils.sysconfig
import distutils.spawn
from distutils.version import LooseVersion
from string import Template

import SCons

import loos_build_config

import scons_support
import platform


EnsureSConsVersion(2, 0)

# ----------------------------------------------------------------------------------------------
# Principal options...
opts = Variables('custom.py')
opts.Add('debug', 'Set to 1 to add -DDEBUG to build', 0)
opts.Add('profile', 'Set to 1 to build the code for profiling', 0)
opts.Add('release', 'Set to 1 to configure for release.', 1)
opts.Add('reparse', 'Set to 1 to regenerate parser-related files.', 0)
opts.Add('pyloos', 'Set to 0 to disable building PyLOOS.', 1)
opts.Add('threads', 'Set to 0 to disable using multithreaded libraries and code', 1)
#opts.Add('docs', 'Set to 0 to disable auto-generation of doxygen documentation', 1)

opts.Add(PathVariable('PREFIX', 'Where to install LOOS', '/opt/LOOS', PathVariable.PathAccept))

opts.Add('BOOST', 'Path to BOOST', '')
opts.Add('BOOST_INCLUDE', 'Path to BOOST Includes', '')
opts.Add('BOOST_LIBPATH', 'Path to BOOST Libraries', '')
opts.Add('BOOST_LIBS', 'Boost libraries to link with', '')

opts.Add('ATLAS_LIBPATH', 'Path to ATLAS Libraries', '')
opts.Add('ATLAS_LIBS', 'Atlas libraries to link with', '')

opts.Add('EIGEN', 'Path to eigen3', '')

opts.Add('NETCDF', 'Path to NetCDF', '')
opts.Add('NETCDF_INCLUDE', 'Path to NetCDF include files', '')
opts.Add('NETCDF_LIBPATH', 'Path to NetCDF libraries', '')
opts.Add('NETCDF_LIBS', 'NetCDF Libraries to link with', '')

opts.Add('PYTHON_PATH', 'Path to Python Modules', '')
opts.Add('PYTHON_INC', 'Include path for Python needed by PyLOOS (if not set, uses the same python as scons)', '')

opts.Add('INCLUDE_PATH', 'Add to include paths before any auto-config', '')
opts.Add('LIBRARY_PATH', 'Add to library paths before any auto-config', '')

scons_support.addDeprecatedOptions(opts)

# If we're using conda, we want to pull in the environment.
# Otherwise, we want the environment mostly cleaned out
swigflags = ['-c++', '-python', '-Wall', '-py3', '-threads']
if "CONDA_PREFIX" in os.environ:
    # Conda has a new enough swig to support the doxygen flag, but most
    # distros don't.
    swigflags.append('-doxygen')
    env = Environment(ENV=os.environ,
                      options=opts,
                      toolpath='.',
                      SWIGFLAGS=swigflags,
                      SHLIBPREFIX=""
                      )
    env["CONDA_PREFIX"] = os.environ["CONDA_PREFIX"]
    env.USING_CONDA = True
else:
    env = Environment(ENV={'PATH': os.environ['PATH']},
                      options=opts,
                      toolpath='.',
                      SWIGFLAGS=swigflags,
                      SHLIBPREFIX=""
                      )
    env.USING_CONDA = False

Help(opts.GenerateHelpText(env))

scons_support.checkForDeprecatedOptions(env)

env.Decider('MD5-timestamp')

# Setup script-builder
script_builder = Builder(action=scons_support.script_builder_python)
env.Append(BUILDERS={'Scripts': script_builder})

# Get more info from environment
PREFIX = env['PREFIX']

# Inject paths (if present)
if 'INCLUDE_PATH' in env:
    env.Append(CPPPATH=env['INCLUDE_PATH'].split(':'))

if 'LIBRARY_PATH' in env:
    env.Append(LIBPATH=env['LIBRARY_PATH'].split(':'))

# ----------------------------------------------------------------------------------------------

cleaning = env.GetOption('clean')


# Autoconf
# TODO: need to update this to handle conda-forge staged recipes
if env.USING_CONDA and platform.system() == "Darwin":
    flag = "-Wl,-rpath," + os.path.join(env["CONDA_PREFIX"], "lib")
    env.Append(LINKFLAGS=flag)

scons_support.AutoConfiguration(env)
pyloos = int(env['pyloos'])

if not pyloos:
    print('***Warning***')
    print('PyLOOS will not be built.  The OMG will not be installed.')


# Compile-flags
debug_opts = '-g -Wall -Wextra -fno-inline'
release_opts = '-O3 -DNDEBUG -Wall -Wno-deprecated'
profile_opts = '-O3 -DNDEBUG -Wall -g'

# Setup the general environment...
env.Prepend(CPPPATH=['#', '#src'])
env.Prepend(LIBPATH=['#', '#src'])
env.Append(LEXFLAGS=['-s'])
env.Append(CPPFLAGS=['-pthread'])
env.Append(LIBS=['pthread'])

# Platform specific build options...
if loos_build_config.host_type == 'Darwin':
    release = platform.release().split('.')
    if int(release[0]) >= 13:    # MacOS 10.9 requires this flag for native compiler
        # Hack to get swig to work with latest 10.9
        env.Append(SWIGFLAGS='-DSWIG_NO_EXPORT_ITERATOR_METHODS')
    env.Append(LINKFLAGS=' -llapack')



if not cleaning:
    # Older version of BOOST will require this definition
    # Note: the version of BOOST requiring this flag is just a guess...
    if LooseVersion(loos_build_config.versions['boost']) < LooseVersion('1_58'):
        env.Append(CCFLAGS='-DBOOST_SPIRIT_USE_PHOENIX_V3=1')

    # This is for BOOST 1.44 and boost 1.45..force using Boost Filesystem v3
    if LooseVersion(loos_build_config.versions['boost']) < LooseVersion('1_46') and LooseVersion(loos_build_config.versions['boost']) >= LooseVersion('1_44'):
        env.Append(CCFLAGS='-DBOOST_FILESYSTEM_VERSION=3')


# Determine what kind of build...

release = int(env['release'])
debug = int(env['debug'])
profile = int(env['profile'])
#docsflag = int(env['docs'])

# If debug is requested, make sure there is no optimization...
if (debug > 0):
   release=0

if int(release):
    env.Append(CCFLAGS=release_opts)
else:
    env.Append(CCFLAGS=debug_opts)

if (debug > 0):
    env.Append(CCFLAGS=(" -DDEBUG=%d" % (debug)))

# Profiling is independent of release/debug status...
if int(profile):
    env.Append(CCFLAGS=profile_opts)
    env.Append(LINKFLAGS=profile_opts)


# Build a revision file to include with LOOS so all tools know what version
# of LOOS they were built with...

scons_support.setupRevision(env)

# Export for subsidiary SConscripts

Export('env')

# ---------------------------------------------------------------------------------------------

# Handle SConscripts and build targets

[loos, loos_python] = SConscript('src/SConscript')
loos_scripts = SConscript('SConscript')
Export('loos')


"""
Handle existing documentation

This is to permit source distributions that include pre-built documentation. If
a docs.prebuilt file is found in the top-level directory, then scons will not
include doxygen sources in the dependency tree (i.e. the docs will not be
rebuilt).  They will also not be included in any cleaning targets.

If a tarball of the documentation is found, then this will be untar'd in place.
It should include a top-level docs.prebuilt file to avoid untar'ing every time.
Tarballs compressed with gzip and bzip2 are recognized, as well as uncompressed
tarballs.

If docs.prebuilt does NOT exist and no tarball is found, then scons will
automatically generate the documentation for most builds, and it will be
included in cleaning.  In addition, install will generate the documentation
"""

"""
if os.path.exists('docs.prebuilt'):
    existing_docs = True
    print('Warning- existing documentation found and will NOT be rebuilt (or cleaned)!')
    print('         Remove docs.prebuilt file to force rebuilding documentation.')
    docs = ['Docs/html/index.html']


else:
    doc_tarballs = glob.glob('loos-*-docs.tar*')
    if doc_tarballs:
        filename = doc_tarballs[0]
        name, extension = os.path.splitext(doc_tarballs[0])
        if extension == '.gz':
            modifier = 'z'
        elif extension == '.bz2':
            modifier = 'j'
        elif extension != '.tar':
            print('Error- unknown compression extension for ', doc_tarballs[0])
            sys.exit(-1)

        print('Warning- existing documentation tarball found.  To force rebuilding of')
        print('         of documentation, remove the tarball and the docs.prebuilt file.')

        if not cleaning:
            print('Unpacking documentation...')
            fnull = open(os.devnull, 'w')
            subprocess.call(['tar', modifier + 'xvf', filename], stdout=fnull)

        existing_docs = True
        docs = ['Docs/html/index.html']

    else:
        existing_docs = False
        docs = env.Doxygen('Doxyfile')
"""

loos_tools = SConscript('Tools/SConscript')

loos_core = loos + loos_scripts


# Automatically setup build targets based on package_list
# Note: excludes Python PyLOOS was not built...

loos_packages = []
for name in loos_build_config.package_list:
    if name == 'Python' and not pyloos:
        continue
    pkg_sc = SConscript('Packages/' + loos_build_config.package_list[name] +
                        '/SConscript')
    env.Alias(name, pkg_sc)
    loos_packages = loos_packages + pkg_sc


# Always install documentation.  Note: html version is hard-coded
"""
env.Command(PREFIX + '/docs/index.html', 'Docs/html/index.html', [
      Delete(PREFIX + '/docs'),
      Copy(PREFIX + '/docs', 'Docs/html'),
      ])
env.AlwaysBuild(PREFIX + '/docs/index.html')
"""


all = loos_tools + loos_scripts + loos_packages

if int(env['pyloos']):
    loos_core = loos_core + loos_python
    all = all + loos_python

loos_tools += loos_core

env.Alias('tools', loos_tools)
env.Alias('core', loos_core)
env.Alias('all', all)
env.Alias('install', PREFIX)

# To get a real "distclean", run "scons -c" and then "scons -c config"
env.Clean('config',
          [
              ".sconsign.dblite",
              ".sconf_temp",
              "config.log"
              ])

env.Default('all')
