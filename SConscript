# SConscript
# (c) 2008 Tod D. Romo
#
# Grossfield Lab
# Department fo Biochemistry & Biophysics
# University of Rochester Medical School
#
#

import sys

Import('env')


library_files = Split('dcd.cpp utils.cpp dcd_utils.cpp AtomicGroup.cpp pdb_remarks.cpp pdb.cpp psf.cpp KernelValue.cpp grammar.yy scanner.ll')
loos = env.Library('loos', library_files)


docs = env.Doxygen('Doxyfile')
examples = SConscript('Examples/SConscript')
tests = SConscript('Tests/SConscript')
tools = SConscript('Tools/SConscript')

env.Alias('docs', docs)
env.Alias('examples', examples)
env.Alias('tests', tests)
env.Alias('tools', tools)
env.Alias('all', loos + examples + tests + docs + tools)

env.Default(loos)
