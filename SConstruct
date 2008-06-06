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

# Principal options...
clos = Options('custom.py')
clos.AddOptions(
	('regenerate', 'Set to 1 to regenerate test outputs', 0),
	('debug', 'Set to 1 to add -DDEBUG to build', 0),
	('release', 'Set to 1 to configure for release.', 0),
	('reparse', 'Set to 1 to regenerate parser-related files.', 0),
)

env = Environment(options = clos, tools = ["default", "doxygen"], toolpath = '.')
Help(clos.GenerateHelpText(env))

# vestigial...
regenerate = env['regenerate']
env['REGENERATE'] = regenerate

reparse = env['reparse']

platform = sys.platform
env['platform'] = platform

# Some rudimentary autoconfish stuff...
if not env.GetOption('clean'):
   conf = Configure(env)
   if not conf.CheckLib('boost_regex'):
      print "***ERROR*** You must have the Boost regular expression libraries installed"
      Exit(1)

   if platform == 'linux2':
      prior = env.get('LIBPATH')
      for dir in ['', '/usr/lib64/atlas', '/usr/lib/atlas', '/usr/local/atlas']:
      	  checks = 1
	  missing = []
	  if dir != "":
	     env.Replace(LIBPATH = [dir])
	     print "Checking for libraries in %s..." % dir
          else:
	     print "Checking for libraries..."

	  if not conf.CheckLib('lapack'):
	     checks = 0
	     missing += ['lapack']
	  if not conf.CheckLib('atlas'):
	     checks = 0
	     missing += ['atlas']
	  if checks:
	     break
      if not checks:
          print "***ERROR*** Missing libraries: ", missing
	  Exit(1)

      if prior != None:
          env['LIBPATH'] = prior + env['LIBPATH']

   env = conf.Finish()

### Compile-flags

debug_opts='-g -Wall -fno-inline'
release_opts='-O3 -DNDEBUG'

# Setup the general environment...
env.Append(CPPPATH = ['#'])
env.Append(LIBPATH = ['#'])
env.Append(LIBS = ['loos', 'boost_regex'])


# Platform specific build options...
if platform == 'darwin':
   env.Append(LINKFLAGS = ' -framework vecLib')

elif platform == 'linux2':
   env.Append(LIBS = ['lapack', 'atlas'])


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



# Export for subsidiary SConscripts
Export('env')

###################################

# Build the LOOS library...
#library_files = Split('dcd.cpp utils.cpp dcd_utils.cpp AtomicGroup.cpp pdb_remarks.cpp pdb.cpp psf.cpp KernelValue.cpp grammar.yy scanner.ll ensembles.cpp')
library_files = Split('dcd.cpp utils.cpp dcd_utils.cpp AtomicGroup.cpp pdb_remarks.cpp pdb.cpp psf.cpp KernelValue.cpp ensembles.cpp dcdwriter.cpp Fmt.cpp')


if int(reparse):
   library_files += ['grammar.yy', 'scanner.ll']
else:
   library_files += ['scanner.cc', 'grammar.cc']


loos = env.Library('loos', library_files)

env.Default(loos)


docs = env.Doxygen('Doxyfile')
examples = SConscript('Examples/SConscript')
tests = SConscript('Tests/SConscript')
tools = SConscript('Tools/SConscript')


# build targets...

env.Alias('docs', docs)
env.Alias('examples', examples)
env.Alias('tests', tests)
env.Alias('tools', tools)

env.Alias('all', loos + examples + tools)
env.Alias('caboodle', loos + examples + tools + tests + docs)

if int(regenerate):
   env.Default('caboodle')
