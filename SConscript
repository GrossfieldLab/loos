# SConscript
# (c) 2008 Tod D. Romo
#
# Grossfield Lab
# Department fo Biochemistry & Biophysics
# University of Rochester Medical School
#
#



Import('env')
cpppath = ['.']

library_files = Split('dcd.cpp utils.cpp dcd_utils.cpp AtomicGroup.cpp pdb_remarks.cpp pdb.cpp psf.cpp KernelValue.cpp grammar.yy scanner.ll')
loos = env.Library('loos', library_files, CPPPATH=cpppath)


docs = env.Doxygen('Doxyfile')
examples = SConscript('Examples/SConscript')
tests = SConscript('Tests/SConscript')

env.Alias('docs', docs)
env.Alias('examples', examples)
env.Alias('tests', tests)
env.Alias('all', loos + examples + tests + docs)

env.Default(loos)
