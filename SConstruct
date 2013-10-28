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
import re
from subprocess import *
from time import strftime
import shutil
import distutils.sysconfig
import distutils.spawn

# Set default path depending on platform...
# Note: This can be reset in custom.py
default_lib_path = '/usr/lib64'

linux_type = 'nonlinux'
platform = sys.platform
if platform == 'linux2':
   # Determine linux variant...
   fv = open('/proc/version', 'r')
   linux_type = fv.read()

   if (re.search("[Uu]buntu", linux_type)):
      default_lib_path = '/usr/lib'


# This is the version-tag for LOOS output
loos_version = '2.1.0'


# Principal options...
clos = Variables('custom.py')
clos.Add('regenerate', 'Set to 1 to regenerate test outputs', 0)
clos.Add('debug', 'Set to 1 to add -DDEBUG to build', 0)
clos.Add('profile', 'Set to 1 to build the code for profiling', 0)
clos.Add('release', 'Set to 1 to configure for release.', 1)
clos.Add('reparse', 'Set to 1 to regenerate parser-related files.', 0)
clos.Add('shared', 'Set to 1 to build a shared LOOS library.', 1)
clos.Add('pyloos', 'Set to 1 to build the python interface to LOOS (requires SWIG).', 0)


clos.Add(PathVariable('LAPACK', 'Path to LAPACK', default_lib_path, PathVariable.PathAccept))
clos.Add(PathVariable('ATLAS', 'Path to ATLAS', default_lib_path + '/atlas', PathVariable.PathAccept))
clos.Add(PathVariable('ATLASINC', 'Path to ATLAS includes', '/usr/include/atlas', PathVariable.PathAccept))
clos.Add(PathVariable('BOOSTLIB', 'Path to BOOST libraries', '', PathVariable.PathAccept))
clos.Add(PathVariable('BOOSTINC', 'Path to BOOST includes', '', PathVariable.PathAccept))
clos.Add('BOOSTSUFFIX', 'Boost Library Name Suffix', '')
clos.Add('CXX', 'C++ Compiler', 'g++')
clos.Add(PathVariable('LIBXTRA', 'Path to additional libraries', '', PathVariable.PathAccept))
clos.Add(PathVariable('PREFIX', 'Path to install LOOS as', '/opt',
                    PathVariable.PathAccept))
clos.Add(PathVariable('NETCDFINC', 'Path to netcdf include files', '', PathVariable.PathAccept))
clos.Add(PathVariable('NETCDFLIB', 'Path to netcdf library files', '', PathVariable.PathAccept))
clos.Add(PathVariable('ALTPATH', 'Additional path to commands', '', PathVariable.PathAccept))
clos.Add(PathVariable('LIBS_OVERRIDE', 'Override linked libs', '', PathVariable.PathAccept))
clos.Add(PathVariable('LIBS_PATHS_OVERRIDE', 'Override paths to libs', '', PathVariable.PathAccept))

# This is a developer setting...  Do not set unless you know what you
# are doing...
clos.Add('REVISION', 'Add build information', loos_version)



env = Environment(options = clos, tools = ["default", "doxygen"], toolpath = '.',SWIGFLAGS=['-c++', '-python', '-Wall'],CPPPATH=[distutils.sysconfig.get_python_inc()],SHLIBPREFIX="")
Help(clos.GenerateHelpText(env))

env.Decider('MD5-timestamp')

# vestigial...
regenerate = env['regenerate']
env['REGENERATE'] = regenerate

reparse = env['reparse']

# export platform to environment...
env['platform'] = platform

env['linux_type'] = linux_type

LAPACK = env['LAPACK']
ATLAS = env['ATLAS']
ATLASINC = env['ATLASINC']
BOOSTLIB = env['BOOSTLIB']
BOOSTINC = env['BOOSTINC']
BOOSTSUFFIX = env['BOOSTSUFFIX']
LIBXTRA = env['LIBXTRA']
PREFIX = env['PREFIX']
ALTPATH = env['ALTPATH']
LIBS_OVERRIDE = env['LIBS_OVERRIDE']
LIBS_PATHS_OVERRIDE = env['LIBS_PATHS_OVERRIDE']
NETCDFINC = env['NETCDFINC']
NETCDFLIB = env['NETCDFLIB']


if ALTPATH != '':
   buildenv = env['ENV']
   path = buildenv['PATH']
   path = ALTPATH + ':' + path
   buildenv['PATH'] = path


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

script_builder = Builder(action = script_builder_python)
env.Append(BUILDERS = {'Scripts' : script_builder})


### Check for swig
swig_location = distutils.spawn.find_executable('swig', env['ENV']['PATH'])



### Autoconf
conf = Configure(env)
if not conf.CheckType('ulong','#include <sys/types.h>\n'):
   conf.env.Append(CCFLAGS = '-DREQUIRES_ULONG')
if not conf.CheckType('uint','#include <sys/types.h>\n'):
   conf.env.Append(CCFLAGS = '-DREQUIRES_UINT')
if conf.CheckLibWithHeader('netcdf', 'netcdf.h', 'c'):    # Should we check C or C++?
   has_netcdf = 1
else:
   has_netcdf = 0
env = conf.Finish()

if (NETCDFINC != '' or NETCDFLIB != ''):
   has_netcdf = 1

env['HAS_NETCDF'] = has_netcdf


### Compile-flags

debug_opts='-g -Wall -Wextra -fno-inline'
release_opts='-O3 -DNDEBUG -Wall'
profile_opts='-pg'

# Setup the general environment...
env.Append(CPPPATH = ['#'])

# Ideally, what's below should be added to the CPPPATH above, but
# doing so causes SCons to scan headers from that directory generating
# implicit dependencies.  SCons seems to mangle these so changing one
# file ends up forcing a complete rebuild.  Setting the include dirs
# directly solves this problem, but it does mean that changes to the
# include files in BOOST and ATLAS will not be picked up by SCons...
if BOOSTINC != '':
   env.Append(CPPFLAGS = ['-I' + BOOSTINC])
env.Append(LIBPATH = ['#', BOOSTLIB, LIBXTRA])
env.Append(LEXFLAGS=['-s'])

LIBS_LINKED_TO = ''
LIBS_PATHS_TO = ''
if (has_netcdf):
   LIBS_LINKED_TO = LIBS_LINKED_TO + ' netcdf'
   env.Append(CPPFLAGS = ['-DHAS_NETCDF'])
   if (NETCDFINC != ''):
      env.Append(CPPFLAGS = ['-I' + NETCDFINC])
   if (NETCDFLIB != ''):
      env.Append(LIBPATH = [NETCDFLIB])


# Platform specific build options...
if platform == 'darwin':
   env.Append(LINKFLAGS = ' -framework Accelerate')
elif platform == 'freebsd8':
   LIBS_LINKED_TO = LIBS_LINKED_TO + ' lapack blas'
elif platform == 'linux2':
   noatlas = 0

   ### Note for OpenSUSE and Ubuntu...
   ### Older versions of those distros may require the gfortran
   ### package be linked in.  If you see strange link errors for
   ### unresolved symbols, try adding "gfortran" to the LIBS list
   ### for your OS below...

   # OpenSUSE doesn't have an atlas package, so use native lapack/blas
   if (re.search("[Ss][Uu][Ss][Ee]", linux_type)):
      LIBS_LINKED_TO = LIBS_LINKED_TO + ' lapack blas'

   elif (re.search("[Uu]buntu", linux_type)):
      LIBS_LINKED_TO = LIBS_LINKED_TO + ' lapack_atlas lapack atlas blas'
      LIBS_PATHS_TO = ATLAS + ' ' + LAPACK

   elif (re.search("[Dd]ebian", linux_type)):
      LIBS_LINKED_TO = LIBS_LINKED_TO + ' atlas lapack blas'
      LIBS_PATHS_TO = ATLAS + ' ' + LAPACK

   else:
      LIBS_LINKED_TO = LIBS_LINKED_TO + ' atlas lapack f77blas'
      LIBS_PATHS_TO = ATLAS + ' ' + LAPACK

   



# CYGWIN does not have an atlas package, so use lapack/blas instead
elif (platform == 'cygwin'):
   LIBS_LINKED_TO = LIBS_LINKED_TO + ' lapack blas'
   LIB_PATHS_TO = 'LAPACK'
   if (BOOSTSUFFIX == ''):
      BOOSTSUFFIX='-mt'

if LIBS_OVERRIDE != '':
   LIBS_LINKED_TO = LIBS_OVERRIDE

if LIBS_PATHS_OVERRIDE != '':
   LIBS_PATHS_TO = LIBS_PATHS_OVERRIDE


# Handle boost after OS-specific options so suffix will be correct
env.Append(LIBS = ['boost_regex' + BOOSTSUFFIX,
                   'boost_program_options' + BOOSTSUFFIX,
                   'boost_thread' + BOOSTSUFFIX,
                   'boost_system' + BOOSTSUFFIX
                   ])


env.Append(LIBS = Split(LIBS_LINKED_TO))
env.Append(LIBPATH = Split(LIBS_PATHS_TO))



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


# Allow overrides from environment...
if os.environ.has_key('CXX'):
   CXX = os.environ['CXX']
   print "Changing default compiler to ", CXX
   env['CXX'] = CXX

if os.environ.has_key('CCFLAGS'):
   CCFLAGS = os.environ['CCFLAGS']
   print "Changing CCFLAGS to ", CCFLAGS
   env['CCFLAGS'] = CCFLAGS

# Divine the current revision...
revision = ''
if env['REVISION'] == '':
   revision = Popen(["svnversion"], stdout=PIPE).communicate()[0]
   revision = revision.rstrip("\n")
else:
   revision = env['REVISION']
      
revision = revision + " " + strftime("%y%m%d")

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

# Export for subsidiary SConscripts

Export('env')



###################################

[loos, loos_python, loos_scripts] = SConscript('SConscript')
Export('loos')

docs = env.Doxygen('Doxyfile')
tests = SConscript('Tests/SConscript')
tools = SConscript('Tools/SConscript')
elastic_networks_package = SConscript('Packages/ElasticNetworks/SConscript')
h_tools = SConscript('Packages/HydrogenBonds/SConscript')
#g_tools = SConscript('Packages/DensityTools/SConscript')
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

env.Alias('lib', loos + loos_scripts)
env.Alias('docs', docs)
env.Alias('tests', tests)
env.Alias('tools', tools)
env.Alias('pyloos_only', loos_python + loos_scripts)

all_target = loos + tools + all_packages + loos_scripts
if env['pyloos'] == '1':
   if swig_location == None:
      print 'Error- Swig (v2.0+) is required to build PyLOOS'
      sys.exit(-1)
   all_target = all_target + loos_python

env.Alias('all', all_target)
env.Alias('caboodle', loos + tools + all_packages + tests + docs + loos_scripts + loos_python)
env.Alias('user', user_package)


env.Alias('install', PREFIX)

if int(regenerate):
   env.Default('caboodle')
else:
   env.Default('all')
