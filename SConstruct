#!/usr/bin/env python
#  This file is part of LOOS.
#
#  LOOS (Lightweight Object-Oriented Structure library)
#  Copyright (c) 2008, Tod D. Romo
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
	('release', 'Set to 1 to configure for release.', 0),
	('reparse', 'Set to 1 to regenerate parser-related files.', 0),
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



### Compile-flags

debug_opts='-g -Wall -Wextra -fno-inline'
release_opts='-O3 -DNDEBUG -Wall'

# Setup the general environment...
env.Append(CPPPATH = ['#', BOOSTINC])
env.Append(LIBPATH = ['#', BOOSTLIB, LIBXTRA])
env.Append(LIBS = ['loos', BOOSTREGEX, BOOSTPO])
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
release = env['release']
if int(release):
    env.Append(CCFLAGS=release_opts)
else:
    env.Append(CCFLAGS=debug_opts)

debug = env['debug']
if int(debug):
   if int(release):
      print "***ERROR*** You cannot have a release with debugging code included."
      Exit(1)
   env.Append(CCFLAGS=" -DDEBUG=$debug")


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

loos = SConscript('SConscript')



docs = env.Doxygen('Doxyfile')
tests = SConscript('Tests/SConscript')
tools = SConscript('Tools/SConscript')


# build targets...

env.Alias('lib', loos)
env.Alias('docs', docs)
env.Alias('tests', tests)
env.Alias('tools', tools)

env.Alias('all', loos + tools)
env.Alias('caboodle', loos + tools + tests + docs)



if int(regenerate):
   env.Default('caboodle')
else:
   env.Default('all')
