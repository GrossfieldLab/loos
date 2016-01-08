#!/usr/bin/env python
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
import glob
import platform
import re
from subprocess import *
from time import strftime
import shutil
import distutils.sysconfig
import distutils.spawn
from string import Template

import SCons
import loos_build_config

from scons_support import *


EnsureSConsVersion(2,0)

# ----------------------------------------------------------------------------------------------


# Principal options...
opts = Variables('custom.py')
opts.Add('debug', 'Set to 1 to add -DDEBUG to build', 0)
opts.Add('profile', 'Set to 1 to build the code for profiling', 0)
opts.Add('release', 'Set to 1 to configure for release.', 1)
opts.Add('reparse', 'Set to 1 to regenerate parser-related files.', 0)
opts.Add('pyloos', 'Set to 0 to disable building PyLOOS.', 1)
opts.Add('threads', 'Set to 0 to disable using multithreaded libraries and code', 1)

opts.Add(PathVariable('PREFIX', 'Where to install LOOS', '/opt/LOOS', PathVariable.PathAccept))

opts.Add('BOOST', 'Path to BOOST', '')
opts.Add('BOOST_INCLUDE', 'Path to BOOST Includes', '')
opts.Add('BOOST_LIBPATH', 'Path to BOOST Libraries', '')
opts.Add('BOOST_LIBS', 'Boost libraries to link with', '')

opts.Add('ATLAS_LIBPATH', 'Path to ATLAS Libraries', '')
opts.Add('ATLAS_LIBS', 'Atlas libraries to link with', '')

opts.Add('NETCDF', 'Path to NetCDF', '')
opts.Add('NETCDF_INCLUDE', 'Path to NetCDF include files', '')
opts.Add('NETCDF_LIBPATH', 'Path to NetCDF libraries', '')
opts.Add('NETCDF_LIBS', 'NetCDF Libraries to link with', '')

opts.Add('PYTHON_PATH', 'Path to Python Modules', '')

addDeprecatedOptions(opts)

### Uncomment this version to have a semi-clean build environment 
#env = Environment(ENV = {'PATH' : os.environ['PATH']}, options = opts, tools = ["default", "doxygen"], toolpath = '.', SWIGFLAGS=['-c++', '-python', '-Wall'],SHLIBPREFIX="")

### Uncomment this line to bring the full user environment into the build environment
env = Environment(ENV = os.environ, options = opts, tools = ["default", "doxygen"], toolpath = '.', SWIGFLAGS=['-c++', '-python', '-Wall'],SHLIBPREFIX="")




Help(opts.GenerateHelpText(env))

checkForDeprecatedOptions(env)

env.Decider('MD5-timestamp')

# Setup script-builder
script_builder = Builder(action = script_builder_python)
env.Append(BUILDERS = {'Scripts' : script_builder})



### Get more info from environment
PREFIX = env['PREFIX']

# ----------------------------------------------------------------------------------------------

### Autoconf
AutoConfiguration(env)
pyloos = int(env['pyloos'])

if not pyloos:
    print '***Warning***'
    print 'PyLOOS will not be built.  The OMG will not be installed.'


### Compile-flags

debug_opts='-g -Wall -Wextra -fno-inline'
release_opts='-O3 -DNDEBUG -Wall'
profile_opts='-O3 -DNDEBUG -Wall -g'

# Setup the general environment...
env.Prepend(CPPPATH = ['#'])
env.Prepend(LIBPATH = ['#'])
env.Append(LEXFLAGS=['-s'])

# Platform specific build options...
if loos_build_config.host_type == 'Darwin':
    release = platform.release().split('.')
    if int(release[0]) >= 13:    # MacOS 10.9 requires this flag for native compiler
        env.Append(CCFLAGS = '--std=c++0x -Wno-deprecated-register -D__ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES=0')
        # Hack to get swig to work with latest 10.9
        env.Append(SWIGFLAGS = '-DSWIG_NO_EXPORT_ITERATOR_METHODS')
    env.Append(LINKFLAGS = ' -framework Accelerate')

# Determine what kind of build...

release = int(env['release'])
debug = int(env['debug'])
profile = int(env['profile'])

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

setupRevision(env)

# Export for subsidiary SConscripts

Export('env')

# ---------------------------------------------------------------------------------------------

### Handle SConscripts and build targets

[loos, loos_python, loos_scripts] = SConscript('SConscript')
Export('loos')

docs = env.Doxygen('Doxyfile')
loos_tools = SConscript('Tools/SConscript')


# Automatically setup build targets based on package_list
# Note: excludes Python PyLOOS was not built...

loos_packages = []
for name in loos_build_config.package_list:
    if name == 'Python' and not pyloos:
        continue
    pkg_sc = SConscript('Packages/' + loos_build_config.package_list[name] + '/SConscript')
    env.Alias(name, pkg_sc)
    loos_packages = loos_packages + pkg_sc


### Special handling for pre-packaged documentation...

env.Command(PREFIX + '/docs/main.html', [], [
      Delete(PREFIX + '/docs'),
      Copy(PREFIX + '/docs', 'Docs'),
      ])
env.AlwaysBuild(PREFIX + '/docs/main.html')


loos_core = loos + loos_scripts
all = loos_tools + loos_scripts + loos_packages

if int(env['pyloos']):
    loos_core = loos_core + loos_python
    all = all + loos_python

loos_tools += loos_core
    
env.Alias('tools', loos_tools)
env.Alias('core', loos_core)
env.Alias('docs', docs)
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
