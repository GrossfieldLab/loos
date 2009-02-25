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

# This can be reset in custom.py
default_lib_path = '/usr/lib64'

# Principal options...
clos = Options('custom.py')
clos.AddOptions(
	('regenerate', 'Set to 1 to regenerate test outputs', 0),
	('debug', 'Set to 1 to add -DDEBUG to build', 0),
        ('profile', 'Set to 1 to build the code for profiling', 0),
	('release', 'Set to 1 to configure for release.', 0),
	('reparse', 'Set to 1 to regenerate parser-related files.', 0),
        ('shared', 'Set to 1 to build a shared LOOS library.', 0)
)

clos.Add(PathOption('LAPACK', 'Path to LAPACK', '', PathOption.PathAccept))
clos.Add(PathOption('ATLAS', 'Path to ATLAS', default_lib_path + '/atlas', PathOption.PathAccept))
clos.Add(PathOption('ATLASINC', 'Path to ATLAS includes', '/usr/include/atlas', PathOption.PathAccept))
clos.Add(PathOption('BOOSTLIB', 'Path to BOOST libraries', '', PathOption.PathAccept))
clos.Add(PathOption('BOOSTINC', 'Path to BOOST includes', '', PathOption.PathAccept))
clos.Add('BOOSTREGEX', 'Boost regex library name', 'boost_regex', PathOption.PathAccept)
clos.Add('BOOSTPO', 'Boost program options library name', 'boost_program_options', PathOption.PathAccept)
clos.Add('CXX', 'C++ Compiler', 'g++')
clos.Add(PathOption('LIBXTRA', 'Path to additional libraries', '', PathOption.PathAccept))
clos.Add(PathOption('PREFIX', 'Path to install LOOS as', '/opt/loos', PathOption.PathAccept))



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
env.Append(CPPPATH = ['#', BOOSTINC])
env.Append(LIBPATH = ['#', BOOSTLIB, LIBXTRA])
env.Append(LIBS = [BOOSTREGEX, BOOSTPO])
env.Append(LEXFLAGS=['-s'])

# Platform specific build options...
if platform == 'darwin':
   env.Append(LINKFLAGS = ' -framework vecLib')
else:
   if platform == 'linux2':
      env.Append(LIBS = ['lapack', 'atlas'])
      env.Append(LIBPATH = [LAPACK, ATLAS])
      env.Append(CPPPATH = [ATLASINC]) 


# Determine what kind of build...
# No option implies debugging, but only an explicit debug defines
# the DEBUG symbol...  Yes, it's a bit obtuse, but it allows
# you to control the level of debugging output through the
# DEBUG definition...

release = env['release']
debug = env['debug']
profile = env['profile']

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

# Export for subsidiary SConscripts
Export('env')

###################################


###################################

if int(env['shared']):
   env['LD_LIBRARY_PATH'] = "."

loos = SConscript('SConscript')



docs = env.Doxygen('Doxyfile')
tests = SConscript('Tests/SConscript')
tools = SConscript('Tools/SConscript')
nm_tools = SConscript('Tools/ElasticNetworks/SConscript')

# Special handling for docs installation...
docs_inst = env.InstallAs(PREFIX + '/share/loos', 'Docs')
Depends(docs_inst, 'Docs/index.html')


# build targets...

env.Alias('lib', loos)
env.Alias('docs', docs)
env.Alias('tests', tests)
env.Alias('tools', tools + nm_tools)

env.Alias('all', loos + tools + nm_tools)
env.Alias('caboodle', loos + tools + nm_tools + tests + docs)

env.Alias('install', ['lib_install', 'tools_install', 'nm_tools_install', docs_inst] )

if int(regenerate):
   env.Default('caboodle')
else:
   env.Default('all')
