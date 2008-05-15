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

library_files = Split('dcd.cpp utils.cpp dcd_utils.cpp AtomicGroup.cpp pdb_remarks.cpp pdb.cpp psf.cpp KernelValue.cpp Parser.cpp Tokenizer.cpp')
loos = env.Library('loos', library_files, CPPPATH=cpppath)


SConscript(['Examples/SConscript', 'Tests/SConscript'])
