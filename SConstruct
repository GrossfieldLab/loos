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

clos.Add(PathOption('LAPACK', 'Path to LAPACK', default_lib_path, PathOption.PathAccept))
clos.Add(PathOption('ATLAS', 'Path to ATLAS', default_lib_path + '/atlas', PathOption.PathAccept))
clos.Add(PathOption('BOOSTLIB', 'Path to BOOST libraries', default_lib_path, PathOption.PathAccept))
clos.Add(PathOption('BOOSTINC', 'Path to BOOST includes', '/usr/include', PathOption.PathAccept))
clos.Add('BOOSTREGEX', 'Boost regex library name', 'boost_regex', PathOption.PathAccept)

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
BOOSTLIB = env['BOOSTLIB']
BOOSTINC = env['BOOSTINC']
BOOSTREGEX = env['BOOSTREGEX']



### Compile-flags

debug_opts='-g -Wall -fno-inline'
release_opts='-O3 -DNDEBUG -Wall'

# Setup the general environment...
env.Append(CPPPATH = ['#', BOOSTINC])
env.Append(LIBPATH = ['#'])
env.Append(LIBS = ['loos', BOOSTREGEX])
env.Append(LEXFLAGS=['-s'])

# Special handling of lib-paths to get only unique paths...
libpaths = { }
libpaths[BOOSTLIB] = 1


# Platform specific build options...
if platform == 'darwin':
   env.Append(LINKFLAGS = ' -framework vecLib')
else:
   if platform == 'linux2':
      libpaths[ATLAS] = 1
      libpaths[LAPACK] = 1
      env.Append(LIBS = ['lapack', 'atlas'])

# Now pull out unique lib paths...
uniquelibs = libpaths.keys()
env.Append(LIBPATH =[uniquelibs])


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

# Build the LOOS library...
library_files = Split('dcd.cpp utils.cpp dcd_utils.cpp pdb_remarks.cpp pdb.cpp psf.cpp KernelValue.cpp ensembles.cpp dcdwriter.cpp Fmt.cpp')
library_files += Split('AtomicGroup.cpp AG_numerical.cpp AG_linalg.cpp Geometry.cpp amber.cpp amber_traj.cpp sfactories.cpp')


if int(reparse):
   library_files += ['grammar.yy', 'scanner.ll']
else:
   library_files += ['scanner.cc', 'grammar.cc']


loos = env.Library('loos', library_files)

env.Default(loos)


docs = env.Doxygen('Doxyfile')
tests = SConscript('Tests/SConscript')
tools = SConscript('Tools/SConscript')


# build targets...

env.Alias('docs', docs)
env.Alias('tests', tests)
env.Alias('tools', tools)

env.Alias('all', loos + tools)
env.Alias('caboodle', loos + tools + tests + docs)

if int(regenerate):
   env.Default('caboodle')
