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



# This can be reset in custom.py
default_lib_path = '/usr/lib64'


# This is the version-tag for LOOS output
loos_version = '1.5.4'


# Principal options...
clos = Variables('custom.py')
clos.Add('regenerate', 'Set to 1 to regenerate test outputs', 0)
clos.Add('debug', 'Set to 1 to add -DDEBUG to build', 0)
clos.Add('profile', 'Set to 1 to build the code for profiling', 0)
clos.Add('release', 'Set to 1 to configure for release.', 1)
clos.Add('reparse', 'Set to 1 to regenerate parser-related files.', 0)
clos.Add('shared', 'Set to 1 to build a shared LOOS library.', 0)
clos.Add('tag', 'Set to 1 to add a revision tag in the code.', 1)


clos.Add(PathVariable('LAPACK', 'Path to LAPACK', default_lib_path, PathVariable.PathAccept))
clos.Add(PathVariable('ATLAS', 'Path to ATLAS', default_lib_path + '/atlas', PathVariable.PathAccept))
clos.Add(PathVariable('ATLASINC', 'Path to ATLAS includes', '/usr/include/atlas', PathVariable.PathAccept))
clos.Add(PathVariable('BOOSTLIB', 'Path to BOOST libraries', '', PathVariable.PathAccept))
clos.Add(PathVariable('BOOSTINC', 'Path to BOOST includes', '', PathVariable.PathAccept))
clos.Add('BOOSTREGEX', 'Boost regex library name', 'boost_regex')
clos.Add('BOOSTPO', 'Boost program options library name', 'boost_program_options')
clos.Add('CXX', 'C++ Compiler', 'g++')
clos.Add(PathVariable('LIBXTRA', 'Path to additional libraries', '', PathVariable.PathAccept))
clos.Add(PathVariable('PREFIX', 'Path to install LOOS as', '/opt',
                    PathVariable.PathAccept))

# This is a developer setting...  Do not set unless you know what you
# are doing...
clos.Add('REVISION', 'Add build information', loos_version)



env = Environment(options = clos, tools = ["default", "doxygen"], toolpath = '.')
Help(clos.GenerateHelpText(env))

# vestigial...
regenerate = env['regenerate']
env['REGENERATE'] = regenerate

reparse = env['reparse']

platform = sys.platform
env['platform'] = platform

LAPACK = env['LAPACK']
ATLAS = env['ATLAS']
ATLASINC = env['ATLASINC']
BOOSTLIB = env['BOOSTLIB']
BOOSTINC = env['BOOSTINC']
BOOSTREGEX = env['BOOSTREGEX']
BOOSTPO = env['BOOSTPO']
LIBXTRA = env['LIBXTRA']
PREFIX = env['PREFIX']


### Autoconf
if not env.GetOption('clean'):
   conf = Configure(env)
   if not conf.CheckType('ulong','#include <sys/types.h>\n'):
      conf.env.Append(CCFLAGS = '-DREQUIRES_ULONG')
   if not conf.CheckType('uint','#include <sys/types.h>\n'):
      conf.env.Append(CCFLAGS = '-DREQUIRES_UINT')
   env = conf.Finish()


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
env.Append(LIBS = [BOOSTREGEX, BOOSTPO])
env.Append(LEXFLAGS=['-s'])

# Platform specific build options...
if platform == 'darwin':
   env.Append(LINKFLAGS = ' -framework vecLib')
else:
   if platform == 'linux2':
      noatlas = 0

      # Determine linux variant
      fv = open('/proc/version', 'r')
      f = fv.read()

      # OpenSUSE doesn't have an atlas package, so use native lapack/blas
      if (re.search("[Ss][Uu][Ss][Ee]", f)):
         env.Append(LIBS = ['lapack', 'blas', 'gfortran'])
         env.Append(LIBPATH = [LAPACK])

      # Ubuntu requires gfortran
      elif (re.search("[Uu]buntu", f)):
         env.Append(LIBS = ['lapack', 'atlas', 'gfortran'])
         env.Append(LIBPATH = [LAPACK, ATLAS])

      # Fedora or similar
      else:
         env.Append(LIBS = ['lapack', 'atlas'])
         env.Append(LIBPATH = [LAPACK, ATLAS])

      #env.Append(CPPPATH = [ATLASINC])       # See above...
      if ATLASINC != '':
         env.Append(CPPFLAGS = ['-I' + ATLASINC])

# Determine what kind of build...
# No option implies debugging, but only an explicit debug defines
# the DEBUG symbol...  Yes, it's a bit obtuse, but it allows
# you to control the level of debugging output through the
# DEBUG definition...

release = int(env['release'])
debug = int(env['debug'])
profile = int(env['profile'])
tag = int(env['tag'])

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

if tag:
   # Divine the current revision...
   revision = ''
   if env['REVISION'] == '':
      revision = Popen(["svnversion"], stdout=PIPE).communicate()[0]
      revision = revision.rstrip("\n")
   else:
      revision = env['REVISION']
      
   revision = revision + " " + strftime("%y%m%d")
   revstr = " -DREVISION=\'\"" + revision + "\"\'"
   env.Append(CCFLAGS=revstr)

### Special handling for pre-packaged documentation...
def DocsInstaller(target, source, env):
   # Get the path to the installed docs dir...
   trgname = target[0].rstr()
   docspath = os.path.split(trgname)
   docsdir = docspath[0]

   # Get path to the source docs dir
   srcpath = os.path.split(source[0].rstr())
   srcdir = os.path.join(srcpath[0], 'Docs')

   shutil.rmtree(docsdir)
   shutil.copytree(srcdir, docsdir)

# Installed docs depend on the docs.built file
env.Command(PREFIX + '/docs/main.html', 'docs.built', DocsInstaller)


# Export for subsidiary SConscripts




Export('env')



###################################

if int(env['shared']):
   env['LD_LIBRARY_PATH'] = "."

loos = SConscript('SConscript')



docs = env.Doxygen('Doxyfile')
tests = SConscript('Tests/SConscript')
tools = SConscript('Tools/SConscript')
nm_tools = SConscript('Tools/ElasticNetworks/SConscript')


# build targets...

env.Alias('lib', loos)
env.Alias('docs', docs)
env.Alias('tests', tests)
env.Alias('tools', tools + nm_tools)

env.Alias('all', loos + tools + nm_tools)
env.Alias('caboodle', loos + tools + nm_tools + tests + docs)

env.Alias('install', PREFIX)

if int(regenerate):
   env.Default('caboodle')
else:
   env.Default('all')
