# Top-level SConstruct
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

env = Environment(tools = ["default", "doxygen"], toolpath = '.')
env.Append(CPPPATH = '#')
env.Append(LIBPATH = '#')
env.Append(LIBS = ['loos', 'boost_regex'])

if sys.platform == 'darwin':
   env.Append(LINKFLAGS = ' -framework vecLib')
elif sys.platform == 'linux':
   env.Append(LIBS = ['lapack', 'atlas'])


# Determine what kind of build...
release=ARGUMENTS.get('release', 0)
if int(release):
    env.Append(CCFLAGS=release_opts)
else:
    env.Append(CCFLAGS=debug_opts)

Export('env')


library_files = Split('dcd.cpp utils.cpp dcd_utils.cpp AtomicGroup.cpp pdb_remarks.cpp pdb.cpp psf.cpp KernelValue.cpp grammar.yy scanner.ll')
loos = env.Library('loos', library_files)

env.Default(loos)


docs = env.Doxygen('Doxyfile')
examples = SConscript('Examples/SConscript')
tests = SConscript('Tests/SConscript')

env.Alias('docs', docs)
env.Alias('examples', examples)
env.Alias('tests', tests)

env.Alias('all', loos + examples + tests + docs)

