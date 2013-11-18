#!/usr/bin/env python
#  This file is part of LOOS.
#
#  LOOS (Lightweight Object-Oriented Structure library)
#  Copyright (c) 2013, Tod D. Romo
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


# This file contains support code for SCons for building LOOS

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


# Attempt to canonicalize system name, type and other related info...
def canonicalizeSystem():
    loos_build_config.host_type = platform.system()

    if os.path.isdir('/usr/lib64'):
        loos_build_config.default_lib_path = '/usr/lib64'
    else:
        loos_build_config.default_lib_path = '/usr/lib'

# Detect CYGWIN & canonicalize linux type, setting defaults...
    if (re.search("(?i)cygwin", loos_build_config.host_type)):
        loos_build_config.host_type = 'Cygwin'
        loos_build_config.suffix = 'dll.a'
    elif (loos_build_config.host_type == 'Linux'):
        # Determine linux variant...
        linux_type = platform.platform()
        
        if (re.search("(?i)ubuntu", linux_type)):
            linux_type = 'ubuntu'
        elif (re.search("(?i)suse", linux_type)):
            linux_type = 'suse'
        elif (re.search("(?i)debian", linux_type)):
            linux_type = 'debian'

        loos_build_config.linux_type = linux_type

    # MacOS is special (of course...)
    elif (loos_build_config.host_type == 'Darwin'):
        loos_build_config.suffix = 'dylib'




### Create a revision file for linking against.
def setupRevision(env):

    # Divine the current revision...
    revision = loos_build_config.loos_version + " " + strftime("%y%m%d")

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


### Let environment variables override or modify some build paramaters...
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

   libpaths = env['LIBPATH']
   libpaths.pop(0)

   cpppaths = env['CPPPATH']
   cpppaths.pop(0)

   if not 'install' in SCons.Script.COMMAND_LINE_TARGETS:
       toolpath = '$LOOS/Tools:' + ':'.join(['$LOOS/Packages/' + s for s in [loos_build_config.package_list[i] for i in loos_build_config.package_list]])
       loos_dir = env.Dir('.').abspath
       libpaths.insert(0, loos_dir)
       cpppaths.insert(0, loos_dir)
       loos_pythonpath = loos_dir

   else:
       loos_dir = env['PREFIX']
       toolpath = loos_dir + '/bin'
       libpaths.insert(0, loos_dir + '/lib')
       loos_pythonpath = loos_dir + '/lib'
       

   file = open(str(source[0]), 'r')
   script = file.read()
   script_template = Template(script)
   script = script_template.substitute(loos_path = loos_dir,
                                       tool_path = toolpath,
                                       libpath = ':'.join(libpaths),
                                       cpppath = ':'.join(cpppaths),
                                       linkflags = env['LINKFLAGS'],
                                       libs = ':'.join(env['LIBS']),
                                       ccflags = env['CCFLAGS'],
                                       loos_cxx = env['CXX'],
                                       loos_pythonpath = loos_pythonpath)

   outfile = open(str(target[0]), 'w')
   outfile.write(script)

   return None



# Verify that we have swig and it's v2.0+
# Returns the path to swig
def CheckForSwig(conf):
    conf.Message('Checking for Swig v2.0+ ...')
    swig_location = distutils.spawn.find_executable('swig', conf.env['ENV']['PATH'])
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



# See if a library requires another to link...
def CheckAtlasRequires(conf, name, lib, required):

    conf.Message('Checking if %s requires %s ... ' % (name, required))
    lastLIBS = conf.env['LIBS']
    conf.env.Append(LIBS=lib)

#    test_code = """
#extern "C"{void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*#, double*, int*, int*);}
#int main(int argc, char *argv[]) { char C[1]; double D[1];int I[1];dgesvd_(C, C, I, I, D, I, D, D, I, #D, I, D, I, I); }
#"""
    test_code="int main(){return(0);}"

    result = conf.TryLink(test_code, '.cpp')
    if not result:
        conf.env.Append(LIBS=required)
        result = conf.TryLink(test_code, '.cpp')
        conf.env.Replace(LIBS=lastLIBS)
        if not result:
            conf.Result('fail')
            return()
        conf.Result('yes')
        return(lib, required)

    conf.env.Replace(LIBS=lastLIBS)
    conf.Result('no')
    return(lib)

    

# Check for existince of boost library with various naming variants
# Will return a tuple containing the correct name and a flag indicating
# whether this is the threaded or non-threaded version.
def CheckForBoostLibrary(conf, name, path, suffix):
   conf.Message('Checking for Boost library %s...' % name)
   name = 'boost_' + name

   if (os.path.isfile(os.path.join(path, 'lib%s-mt.%s' % (name , suffix)))):
      conf.Result(name + '-mt')
      return(name + '-mt', 1)

   if (os.path.isfile(os.path.join(path, 'lib%s.%s' % (name , suffix)))):
      conf.Result(name)
      return(name, 0)

   def sortByLength(w1,w2):
      return len(w1)-len(w2)

    # Now check for names lib libboost_regex-gcc43-mt.so ...
   files = glob.glob(os.path.join(path, 'lib%s-*-mt.%s' % (name, suffix)))
   files.sort(cmp=sortByLength)
   if files:
      conf.Result(name + '-mt')
      name = os.path.basename(files[0])[3:-(len(suffix)+1)]
      return(name, 1)

   files = glob.glob(os.path.join(path, 'lib%s-*.%s' % (name, suffix)))
   files.sort(cmp=sortByLength)
   if files:
      conf.Result(name)
      name = os.path.basename(files[0])[3:-(len(suffix)+1)]
      return(name, 0)


   conf.Result('missing')
   return('', -1)

            
# Check for version of Boost includes
def CheckBoostHeaderVersion(conf, min_boost_version):
    source_code = """
#include <boost/version.hpp>
#if (((BOOST_VERSION / 100) % 1000) < $version)
#error LOOS require Boost 1.$version or higher
#endif
int main(int argc, char *argv[]) { return(0); }
"""

    st = Template(source_code)
    test_code = st.substitute(version = min_boost_version)

    conf.Message('Checking Boost version... ')
    result = conf.TryLink(test_code, '.cpp')
    if not result:
        conf.Result('too old (use Boost 1.%d+)' % min_boost_version)
        return(0)

    conf.Result('ok')
    return(1)

# Check for presence of a directory
def CheckDirectory(conf, dirname):

    conf.Message('Checking for directory %s...' % dirname)
    if os.path.isdir(dirname):
        conf.Result('yes')
        return(1)
    conf.Result('no')
    return(0)


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
        boost_libpath = loos_build_config.default_lib_path
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

    # The main SConstruct will need to know what this path is for autoconf testing...
    env['BOOST_LIBPATH'] = boost_libpath



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
        netcdf_libpath = loos_build_config.default_lib_path
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


def AutoConfiguration(env):
    if env.GetOption('clean') or env.GetOption('help'):
        env['HAS_NETCDF'] = 1
    else:
        has_netcdf = 0
        conf = env.Configure(custom_tests = { 'CheckForSwig' : CheckForSwig,
                                               'CheckForBoostLibrary' : CheckForBoostLibrary,
                                               'CheckBoostHeaderVersion' : CheckBoostHeaderVersion,
                                               'CheckDirectory' : CheckDirectory,
                                               'CheckAtlasRequires' : CheckAtlasRequires})
        
    
        # Some distros use /usr/lib, others have /usr/lib64.
        # Check to see what's here and prefer lib64 to lib
        if not conf.CheckDirectory('/usr/lib64'):
            if not conf.CheckDirectory('/usr/lib'):
                print 'Fatal error- cannot find your system library directory'
                env.Exit(1)
            default_lib_path = '/usr/lib'
       
        # Now that we know the default library path, setup Boost, NetCDF, and ATLAS
        # based on the environment or custom.py file
        SetupBoostPaths(env)
        SetupNetCDFPaths(env)

        # Only setup ATLAS if we're not on a Mac...
        if loos_build_config.host_type != 'Darwin':
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
        if conf.env['NETCDF_LIBS']:
            netcdf_libs = env['NETCDF_LIBS']
            env.Append(CCFLAGS=['-DHAS_NETCDF'])
            has_netcdf = 1
        else:
            if conf.CheckLibWithHeader('netcdf', 'netcdf.h', 'c'):    # Should we check C or C++?
                netcdf_libs = 'netcdf'
                env.Append(CCFLAGS=['-DHAS_NETCDF'])
                has_netcdf = 1

        conf.env['HAS_NETCDF'] = has_netcdf


        # --- Swig Autoconf (unless user requested NO PyLOOS)
        if int(env['pyloos']):
            if conf.CheckForSwig():
                conf.env['pyloos'] = 1
            else:
                conf.env['pyloos'] = 0

        # --- Boost Autoconf
        if not conf.CheckBoostHeaderVersion(loos_build_config.min_boost_version):
            env.Exit(1)

        if conf.env['BOOST_LIBS']:
            boost_libs = Split(env['BOOST_LIBS'])
        else:
            boost_threaded = -1
            boost_libs = []
            for libname in ['regex', 'thread', 'system', 'program_options']:
                result = conf.CheckForBoostLibrary(libname, env['BOOST_LIBPATH'], loos_build_config.suffix)
                if not result[0]:
                    print 'Error- missing Boost library %s' % libname
                    env.Exit(1)
                if boost_threaded < 0:
                    boost_threaded = result[1]
                elif boost_threaded and not result[1]:
                    print 'Error- Expected threaded boost libraries, but %s is not threaded.' % libname
                    env.Exit(1)
                elif not boost_threaded and result[1]:
                    print 'Error- Expected non-threaded boost libraries, but %s is threaded.' % libname
                    env.Exit(1)
                boost_libs.append(result[0])
    
            env.Append(LIBS = boost_libs)


        # --- Check for ATLAS/LAPACK and how to build

        # MacOS will use accelerate framework, so skip all of this...
        if loos_build_config.host_type != 'Darwin':
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
                    env.Exit(1)
                    
                if (numerics['atlas']):
                    atlas_libs.append('atlas')

                if not (numerics['lapack'] or numerics['atlas']):
                    # Did we miss atlas above because it requires gfortran?
                    if not numerics['atlas'] and (numerics['f77blas'] and numerics['cblas']):
                        atlas_libs = conf.CheckAtlasRequires('atlas', 'atlas' + atlas_libs, 'gfortran')
                        if not atlas_libs:
                            print 'Error- could not figure out how to build with Atlas/lapack'

                    # In some cases, lapack requires blas to link so the above
                    # check will find blas but not lapack
                    elif numerics['blas']:
                        result = conf.CheckAtlasRequires('lapack', 'lapack', 'blas')
                        if result:
                            atlas_libs.append('lapack')
                        else:
                            print 'Error- you must have either Lapack or Atlas installed'
                            env.Exit(1)
                    else:
                        print 'Error- you must have either Lapack or Atlas installed'
                        env.Exit(1)

                if not atlas_libs:
                    print 'Error- could not figure out how to build with Atlas/Lapack'
                    env.Exit(1)

                # Hack to extend list rather than append a list into a list
                for lib in atlas_libs:
                    env.Append(LIBS=lib)

        environOverride(conf)
        env = conf.Finish()
