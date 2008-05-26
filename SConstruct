# (c) 2008 Tod D. Romo
#
# Grossfield Lab
# Department fo Biochemistry & Biophysics
# University of Rochester Medical School
#
#

import sys

### Compile-flags

debug_opts='-g -Wall'
release_opts='-O3'

# Setup the environment...
env = Environment(tools = ["default", "doxygen"], toolpath = '.')
env.Append(CPPPATH = '#')
env.Append(LIBPATH = '#')
env.Append(LIBS = ['loos', 'boost_regex'])

# Platform specific build options...
if sys.platform == 'darwin':
   env.Append(LINKFLAGS = ' -framework vecLib')
elif sys.platform == 'linux2':
   env.Append(LIBS = ['lapack', 'atlas'])
   env.Replace(LIBPATH = ['#', '/usr/lib64/atlas'])    # Not sure why I have to replace
                                                       # rather than append...


# Determine what kind of build...
release=ARGUMENTS.get('release', 0)
if int(release):
    env.Append(CCFLAGS=release_opts)
else:
    env.Append(CCFLAGS=debug_opts)

# Export for subsidiary SConscripts
Export('env')

# Build the LOOS library...
library_files = Split('dcd.cpp utils.cpp dcd_utils.cpp AtomicGroup.cpp pdb_remarks.cpp pdb.cpp psf.cpp KernelValue.cpp grammar.yy scanner.ll')
loos = env.Library('loos', library_files)

env.Default(loos)


docs = env.Doxygen('Doxyfile')
examples = SConscript('Examples/SConscript')
tests = SConscript('Tests/SConscript')


# build targets...

env.Alias('docs', docs)
env.Alias('examples', examples)
env.Alias('tests', tests)

env.Alias('all', loos + examples + tests + docs)

