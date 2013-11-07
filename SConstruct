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

import SCons

loos_version = '2.1.0'   # Set to null string to use SVN revision

# Set default path depending on platform...
# Note: This can be reset in custom.py
default_lib_path = '/usr/lib64'



def canonicalizeSystem():
    linux_type = 'nonlinux'
    host_type = platform.system()
    suffix = 'so'
# Detect CYGWIN & canonicalize linux type, setting defaults...
    if (re.search("(?i)cygwin", host_type)):
        host_type = 'Cygwin'
    elif (host_type == 'Linux'):
        # Determine linux variant...
        linux_type = platform.platform()
        
        if (re.search("(?i)ubuntu", linux_type)):
            linux_type = 'ubuntu'
            default_lib_path = '/usr/lib'
        elif (re.search("(?i)suse", linux_type)):
            linux_type = 'suse'
        elif (re.search("(?i)debian", linux_type)):
            linux_type = 'debian'
    elif (host_type == 'Darwin'):
        suffix = 'dylib'
    return(host_type, linux_type, suffix)



def setupRevision(env):
    # Divine the current revision...
    revision = loos_version + " " + strftime("%y%m%d")

    # Now, write this out to a cpp file that can be linked in...this avoids having
    # to recompile everything when building on a new date.  We also rely on SCons
    # using the MD5 checksum to detect changes in the file (even though it's always
    # rewritten)
    revfile = open('revision.cpp', 'w')
    revfile.write('#include <string>\n')
    revfile.write('std::string revision_label = "')
    revfile.write(revision)
    revfile.write('";\n')
    revfile.close()


def environOverride(conf):
    # Allow overrides from environment...
    if 'CXX' in os.environ:
        conf.env.Replace(CXX = os.environ['CXX'])
        print '*** Using compiler ' + os.environ['CXX']
    
    if 'CCFLAGS' in os.environ:
        conf.env.Append(CCFLAGS = os.environ['CCFLAGS'])
        print '*** Appending custom build flags: ' + os.environ['CCFLAGS']
        
    if 'LDFLAGS' in os.environ:
        conf.env.Append(LINKFLAGS = os.environ['LDFLAGS'])
        print '*** Appending custom link flag: ' + os.environ['LDFLAGS']

### Builder for setup scripts

# This copies the environment setup script while changing the directory
# that's used for setting up PATH and [DY]LD_LIBRARY_PATH.  If LOOS
# is being built in a directory, the env script will be setup to use
# the built-in-place distribution.  If LOOS is being installed, then
# it will use the installation directory instead.

def script_builder_python(target, source, env):
   first = target[0]
   target_path = first.get_abspath()
   dir_path = os.path.dirname(target_path)

   command = "sed s@PATH_TO_LOOS@" + dir_path + "@ <" + str(source[0]) + " >" + str(first)
#   print command
   os.system(command)
   return None



def CheckForSwig(conf):
    conf.Message('Checking for Swig v2.0+ ...')
    swig_location = distutils.spawn.find_executable('swig', env['ENV']['PATH'])
    if swig_location == None:
        conf.Result('no')
    else:
        swig_check = Popen([swig_location, "-version"], stdout=PIPE).communicate()[0]
        swig_version = swig_check.split('\n')[1].split(' ')[2]
        swig_major = swig_version.split('.')[0]
        if int(swig_major) < 2:
            conf.Result('no')
            swig_location = None
        else:
            conf.Result('yes')
    return(swig_location)



# See if we need gfortran in order to build code with atlas/lapack
# Returns the full list of libraries used...

def CheckAtlasBuild(conf, libs):
   test_code = """
extern "C"{void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);}
int main(int argc, char *argv[]) { char C[1]; double D[1];int I[1];dgesvd_(C, C, I, I, D, I, D, D, I, D, I, D, I, I); }
"""

   conf.Message('Checking if ATLAS/LAPACK needs gfortran...')
   
   lastLIBS = conf.env['LIBS']
   conf.env.Append(LIBS = libs)
   result = conf.TryLink(test_code, '.cpp')
   if not result:
      conf.env.Append(LIBS = option_libs)
      result = conf.TryLink(test_code)
      if not result:
         conf.Result('error')
         return([])
      conf.Result('yes')
      conf.env.Replace(LIBS = lastLIBS)
      return(libs + 'gfortran')
   conf.Result('no')
   conf.env.Replace(LIBS = lastLIBS)
   return(libs)


# Check for existince of boost library with various naming variants
# Will return a tuple containing the correct name and a flag indicating
# whether this is the threaded or non-threaded version.

def CheckForBoostLibrary(conf, name, path, suffix):
   conf.Message('Checking for Boost library %s...' % name)

   if (os.path.isfile(os.path.join(path, 'lib%s-mt.%s' % (name , suffix)))):
      conf.Result('yes')
      return(name + '-mt', 1)

   if (os.path.isfile(os.path.join(path, 'lib%s.%s' % (name , suffix)))):
      conf.Result('yes')
      return(name, 0)

   def sortByLength(w1,w2):
      return len(w1)-len(w2)

    # Now check for names lib libboost_regex-gcc43-mt.so ...
   files = glob.glob(os.path.join(path, 'lib%s-*-mt.%s' % (name, suffix)))
   files.sort(cmp=sortByLength)
   if files:
      conf.Result('yes')
      name = os.path.basename(files[0])[3:-(len(suffix)+1)]
      return(name, 1)

   files = glob.glob(os.path.join(path, 'lib%s-*.%s' % (name, suffix)))
   files.sort(cmp=sortByLength)
   if files:
      conf.Result('yes')
      name = os.path.basename(files[0])[3:-(len(suffix)+1)]
      return(name, 0)


   conf.Result('no')
   return('', -1)

            

# ----------------------------------------------------------------------------------


def SetupBoostPaths(env):

    global BOOST_LIBS
    global boost
    global boost_include
    global boost_libpath

    BOOST=env['BOOST']
    BOOST_INCLUDE=env['BOOST_INCLUDE']
    BOOST_LIBPATH=env['BOOST_LIBPATH']
    BOOST_LIBS = env['BOOST_LIBS']

    if BOOST == '':
        boost = '/usr'
        boost_include = '/usr/include'
        boost_libpath = '/usr/lib64'
    else:
        boost = BOOST
        boost_include = boost + '/include'
        boost_libpath = boost + '/lib'
        
    if BOOST_INCLUDE:
        boost_include = BOOST_INCLUDE
    if BOOST_LIBPATH:
        boost_libpath= BOOST_LIBPATH
       

    env.MergeFlags({ 'LIBPATH': [boost_libpath]})
    env.MergeFlags({ 'CPPPATH' : [boost_include] })


########################3


def SetupNetCDFPaths(env):
    global NETCDF_LIBS
    global netcdf_include
    global netcdf_libpath

    NETCDF=env['NETCDF']
    NETCDF_INCLUDE=env['NETCDF_INCLUDE']
    NETCDF_LIBPATH=env['NETCDF_LIBPATH']
    NETCDF_LIBS = env['NETCDF_LIBS']
    
    if NETCDF == '':
        netcdf = '/usr'
        netcdf_include = '/usr/include'
        netcdf_libpath = '/usr/lib64'
    else:
        netcdf = NETCDF
        netcdf_include = netcdf + '/include'
        netcdf_libpath = netcdf + '/lib'

    if NETCDF_INCLUDE:
        netcdf_include = NETCDF_INCLUDE
    if NETCDF_LIBPATH:
        netcdf_libpath= NETCDF_LIBPATH


    env.MergeFlags({ 'LIBPATH': [netcdf_libpath]})
    env.MergeFlags({ 'CPPPATH' : [netcdf_include] })


# This is the version-tag for LOOS output
loos_version = '2.1.0'


# Principal options...
opts = Variables('custom.py')
opts.Add('debug', 'Set to 1 to add -DDEBUG to build', 0)
opts.Add('profile', 'Set to 1 to build the code for profiling', 0)
opts.Add('release', 'Set to 1 to configure for release.', 1)
opts.Add('reparse', 'Set to 1 to regenerate parser-related files.', 0)
opts.Add('pyloos', 'Set to 0 to disable building PyLOOS.', 1)

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


### Setup paths...

#--- First, Boost
SetupBoostPaths(env)

#--- And now NETCDF
SetupNetCDFPaths(env)

#--- Now, ATLAS

ATLAS_LIBPATH = env['ATLAS_LIBPATH']
ATLAS_LIBS = env['ATLAS_LIBS']
if not ATLAS_LIBPATH:
    atlas_libpath = '/usr/lib64/atlas'
else:
    atlas_libpath = ATLAS_LIBPATH

env.MergeFlags({ 'LIBPATH': [atlas_libpath] })


### Get more info from environment
PREFIX = env['PREFIX']

### Autoconf
# (don't bother when cleaning)
has_netcdf = 0
pyloos = int(env['pyloos'])

if not env.GetOption('clean'):
    conf = Configure(env, custom_tests = { 'CheckForSwig' : CheckForSwig,
                                           'CheckAtlasBuild' : CheckAtlasBuild,
                                           'CheckForBoostLibrary' : CheckForBoostLibrary })

    if not conf.CheckType('ulong','#include <sys/types.h>\n'):
        conf.env.Append(CCFLAGS = '-DREQUIRES_ULONG')
    if not conf.CheckType('uint','#include <sys/types.h>\n'):
        conf.env.Append(CCFLAGS = '-DREQUIRES_UINT')

# --- NetCDF
    env['HAS_NETCDF'] = 0
    if NETCDF_LIBS:
        netcdf_libs = NETCDF_LIBS
        env.Append(CCFLAGS=['-DHAS_NETCDF'])
        env['HAS_NETCDF'] = 1
    else:
        if conf.CheckLibWithHeader('netcdf', 'netcdf.h', 'c'):    # Should we check C or C++?
            netcdf_libs = 'netcdf'
            env.Append(CCFLAGS=['-DHAS_NETCDF'])
            env['HAS_NETCDF'] = 1

# --- SWIG
    if pyloos:
        if conf.CheckForSwig():
            pyloos = 1
        else:
            pyloos = 0

    env['pyloos'] = pyloos

# --- BOOST
    if BOOST_LIBS:
        boost_libs = Split(BOOST_LIBS)
    else:
        boost_threaded = -1
        boost_libs = []
        for libname in ['boost_regex', 'boost_thread', 'boost_system', 'boost_program_options']:
            result = conf.CheckForBoostLibrary(libname, boost_libpath, 'so')
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

# --- ATLAS
    if ATLAS_LIBS:
        atlas_libs = Split(ATLAS_LIBS)
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


    environOverride(conf)

    env = conf.Finish()

    env.Append(LIBS = boost_libs)
    env.Append(LIBS = atlas_libs)
    


### Compile-flags

debug_opts='-g -Wall -Wextra -fno-inline'
release_opts='-O3 -DNDEBUG -Wall'
profile_opts='-pg'

# Setup the general environment...
env.Prepend(CPPPATH = ['#'])
env.Prepend(LIBPATH = ['#'])
env.Append(LEXFLAGS=['-s'])

if pyloos:
    env.Append(CPPPATH = [distutils.sysconfig.get_python_inc()])

# Platform specific build options...
if host_type == 'Darwin':
    release = platform.release().split('.')
    if int(release[0]) >= 13:    # MacOS 10.9 requires this flag for native compiler
        env.Append(CCFLAGS = '--std=c++0x')
    env.Append(LINKFLAGS = ' -framework Accelerate')

# Determine what kind of build...
# No option implies debugging, but only an explicit debug defines
# the DEBUG symbol...  Yes, it's a bit obtuse, but it allows
# you to control the level of debugging output through the
# DEBUG definition...

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


setupRevision(env)

# Export for subsidiary SConscripts

Export('env')



###################################

[loos, loos_python, loos_scripts] = SConscript('SConscript')
Export('loos')

docs = env.Doxygen('Doxyfile')
#tests = SConscript('Tests/SConscript')
tools = SConscript('Tools/SConscript')
elastic_networks_package = SConscript('Packages/ElasticNetworks/SConscript')
h_tools = SConscript('Packages/HydrogenBonds/SConscript')
convergence_package = SConscript('Packages/Convergence/SConscript')
density_package = SConscript('Packages/DensityTools/SConscript')
user_package = SConscript('Packages/User/SConscript')
python_package = SConscript('Packages/PyLOOS/SConscript')


all_packages = elastic_networks_package + h_tools + convergence_package + density_package

### Special handling for pre-packaged documentation...

env.Command(PREFIX + '/docs/main.html', [], [
      Delete(PREFIX + '/docs'),
      Copy(PREFIX + '/docs', 'Docs'),
      ])
env.AlwaysBuild(PREFIX + '/docs/main.html')


# build targets...

loos_core = loos + loos_scripts
loos_tools = tools
if pyloos:
    loos_core = loos_core + loos_python
    loos_tools = loos_tools + loos_python

env.Alias('core', loos_core)
env.Alias('docs', docs)
#env.Alias('tests', tests)
env.Alias('tools', loos_tools)

all_target = loos + tools + all_packages + loos_scripts
if pyloos:
    all_target = all_target + loos_python

env.Alias('all', all_target)
env.Alias('caboodle', loos + tools + all_packages + docs + loos_scripts + loos_python)
env.Alias('user', user_package)


env.Alias('install', PREFIX)

env.Default('all')
