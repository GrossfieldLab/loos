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


# ----------------------------------------------------------------------------------------------


# Principal options...
opts = Variables('custom.py')
opts.Add('debug', 'Set to 1 to add -DDEBUG to build', 0)
opts.Add('profile', 'Set to 1 to build the code for profiling', 0)
opts.Add('release', 'Set to 1 to configure for release.', 1)
opts.Add('reparse', 'Set to 1 to regenerate parser-related files.', 0)
opts.Add('pyloos', 'Set to 0 to disable building PyLOOS.', 1)

opts.Add(PathVariable('PREFIX', 'Where to install LOOS', '/opt/LOOS', PathVariable.PathAccept))
opts.Add('ALTPATH', 'Additional path to commands', '')

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


env = Environment(options = opts, tools = ["default", "doxygen"], toolpath = '.',SWIGFLAGS=['-c++', '-python', '-Wall'],SHLIBPREFIX="")
Help(opts.GenerateHelpText(env))

env.Decider('MD5-timestamp')

# Setup script-builder
script_builder = Builder(action = script_builder_python)
env.Append(BUILDERS = {'Scripts' : script_builder})

# Get system information
(host_type, linux_type, library_suffix) = canonicalizeSystem()
env['host_type'] = host_type
env['linux_type'] = linux_type


# Setup alternate path to tools
if env['ALTPATH']:
   buildenv = env['ENV']
   path = buildenv['PATH']
   path = ALTPATH + ':' + path
   buildenv['PATH'] = path


### Setup paths...

#--- Now, ATLAS



### Get more info from environment
PREFIX = env['PREFIX']

# ----------------------------------------------------------------------------------------------

### Autoconf
# (don't bother when cleaning)
has_netcdf = 0
pyloos = int(env['pyloos'])
env['HAS_NETCDF'] = 0
default_lib_path = '/usr/lib64'

if not (env.GetOption('clean') or env.GetOption('help')):
    conf = Configure(env, custom_tests = { 'CheckForSwig' : CheckForSwig,
                                           'CheckAtlasBuild' : CheckAtlasBuild,
                                           'CheckForBoostLibrary' : CheckForBoostLibrary,
                                           'CheckBoostHeaderVersion' : CheckBoostHeaderVersion,
                                           'CheckDirectory' : CheckDirectory })

    
    # Some distros use /usr/lib, others have /usr/lib64.
    # Check to see what's here and prefer lib64 to lib
    if not conf.CheckDirectory('/usr/lib64'):
       if not conf.CheckDirectory('/usr/lib'):
          print 'Fatal error- cannot find your system library directory'
          Exit(1)
       default_lib_path = '/usr/lib'

    # Now that we know the default library path, setup Boost, NetCDF, and ATLAS
    # based on the environment or custom.py file
    SetupBoostPaths(env)
    SetupNetCDFPaths(env)

    # Only setup ATLAS if we're not on a Mac...
    if host_type != 'Darwin':
       ATLAS_LIBPATH = env['ATLAS_LIBPATH']
       ATLAS_LIBS = env['ATLAS_LIBS']
       if not ATLAS_LIBPATH:
          atlas_libpath = loos_build_config.default_lib_path + '/atlas'
       else:
          atlas_libpath = ATLAS_LIBPATH

       env.MergeFlags({ 'LIBPATH': [atlas_libpath] })


    # Check for standard typedefs...
    if not conf.CheckType('ulong','#include <sys/types.h>\n'):
        conf.env.Append(CCFLAGS = '-DREQUIRES_ULONG')
    if not conf.CheckType('uint','#include <sys/types.h>\n'):
        conf.env.Append(CCFLAGS = '-DREQUIRES_UINT')


# --- NetCDF Autoconf
    has_netcdf = 0
    if env['NETCDF_LIBS']:
        netcdf_libs = env['NETCDF_LIBS']
        env.Append(CCFLAGS=['-DHAS_NETCDF'])
        has_netcdf = 1
    else:
        if conf.CheckLibWithHeader('netcdf', 'netcdf.h', 'c'):    # Should we check C or C++?
            netcdf_libs = 'netcdf'
            env.Append(CCFLAGS=['-DHAS_NETCDF'])
            has_netcdf = 1

    env['HAS_NETCDF'] = has_netcdf


# --- Swig Autoconf (unless user requested NO PyLOOS)
    if pyloos:
        if conf.CheckForSwig():
            pyloos = 1
        else:
            pyloos = 0

    env['pyloos'] = pyloos

# --- Boost Autoconf
    if env['BOOST_LIBS']:
        boost_libs = Split(env['BOOST_LIBS'])
    else:
        boost_threaded = -1
        boost_libs = []
        for libname in ['regex', 'thread', 'system', 'program_options']:
            result = conf.CheckForBoostLibrary(libname, env['BOOST_LIBPATH'], library_suffix)
            if not result[0]:
                print 'Error- missing Boost library %s' % libname
                Exit(1)
            if boost_threaded < 0:
                boost_threaded = result[1]
            elif boost_threaded and not result[1]:
                print 'Error- Expected threaded boost libraries, but %s is not threaded.' % libname
                Exit(1)
            elif not boost_threaded and result[1]:
                print 'Error- Expected non-threaded boost libraries, but %s is threaded.' % libname
                Exit(1)
            boost_libs.append(result[0])
    
    env.Append(LIBS = boost_libs)

    if not conf.CheckBoostHeaderVersion(loos_build_config.min_boost_version):
        Exit(1)

# --- Check for ATLAS/LAPACK and how to build

    # MacOS will use accelerate framework, so skip all of this...
    if host_type != 'Darwin':
        if env['ATLAS_LIBS']:
            atlas_libs = Split(env['ATLAS_LIBS'])
        else:
            numerics = { 'atlas' : 0,
                         'lapack' : 0,
                         'f77blas' : 0,
                         'cblas' : 0,
                         'blas' : 0 }
            
        
            for libname in numerics.keys():
                if conf.CheckLib(libname, autoadd = 0):
                    numerics[libname] = 1

            atlas_libs = []
            if (numerics['lapack']):
                atlas_libs.append('lapack')
            
            if (numerics['f77blas'] and numerics['cblas']):
                atlas_libs.extend(['f77blas', 'cblas'])
            elif (numerics['blas']):
                atlas_libs.append('blas')
            else:
                print 'Error- you must have some kind of blas installed'
                Exit(1)
                    
            if (numerics['atlas']):
                atlas_libs.append('atlas')

            if not numerics['lapack'] and not numerics['atlas']:
                print 'Error- you must have either LAPACK or Atlas installed'
                Exit(1)

            atlas_libs = conf.CheckAtlasBuild(atlas_libs)
            if not atlas_libs:
                print 'Error- could not figure out how to build.'
                Exit(1)

        env.Append(LIBS = atlas_libs)


    environOverride(conf)
    env = conf.Finish()
    


### Compile-flags

debug_opts='-g -Wall -Wextra -fno-inline'
release_opts='-O3 -DNDEBUG -Wall'
profile_opts='-pg'

# Setup the general environment...
env.Prepend(CPPPATH = ['#'])
env.Prepend(LIBPATH = ['#'])
env.Append(LEXFLAGS=['-s'])

# Include Python if building PyLOOS
if pyloos:
    env.Append(CPPPATH = [distutils.sysconfig.get_python_inc()])

# Platform specific build options...
if host_type == 'Darwin':
    release = platform.release().split('.')
    if int(release[0]) >= 13:    # MacOS 10.9 requires this flag for native compiler
        env.Append(CCFLAGS = '--std=c++0x')
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
   env.Append(CCFLAGS=" -DDEBUG=$debug")

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

if pyloos:
    loos_core = loos_core + loos_python
    loos_tools = loos_tools + loos_python

all = loos_tools + loos_scripts + loos_packages

env.Alias('tools', loos_tools)
env.Alias('core', loos_core)
env.Alias('docs', docs)
env.Alias('all', all)
env.Alias('install', PREFIX)

env.Default('all')
