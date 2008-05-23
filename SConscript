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


#dox_target = env.Alias("dox", env.Doxygen("Doxyfile"))
doxy = env.Doxygen('Doxyfile')
print "Doxy[0] = ", doxy[0], "\n"
print "Doxy[1] = ", doxy[1], "\n"

SConscript(['Examples/SConscript', 'Tests/SConscript'])
