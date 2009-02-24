#!/usr/bin/env python


import sys

Import('env')


reparse = env['reparse']
if int(reparse):
   apps = 'grammar.yy scanner.ll '
else:
   apps = 'scanner.cc grammar.cc '


apps = apps + 'dcd.cpp utils.cpp dcd_utils.cpp pdb_remarks.cpp pdb.cpp psf.cpp KernelValue.cpp ensembles.cpp dcdwriter.cpp Fmt.cpp'
apps = apps + ' AtomicGroup.cpp AG_numerical.cpp AG_linalg.cpp Geometry.cpp amber.cpp amber_traj.cpp tinkerxyz.cpp sfactories.cpp'
apps = apps + ' ccpdb.cpp pdbtraj.cpp tinker_arc.cpp'

if int(env['shared']):
   loos = env.SharedLibrary('loos', Split(apps))
else:
   loos = env.Library('loos', Split(apps))


# Handle installation...
PREFIX = env['PREFIX']

# Library(ies)
loos_lib_inst = env.Install(PREFIX + '/lib', loos)

 
# Header files...
hdr = 'loos.hpp'
hdr = 'amber.hpp amber_traj.hpp Atom.hpp AtomicGroup.hpp ccpdb.hpp Coord.hpp'
hdr = hdr + ' cryst.hpp dcd.hpp dcd_utils.hpp dcdwriter.hpp ensembles.hpp Fmt.hpp'
hdr = hdr + ' Geometry.hpp KernelActions.hpp Kernel.hpp KernelStack.hpp'
hdr = hdr + ' KernelValue.hpp loos_defs.hpp loos.hpp LoosLexer.hpp Matrix44.hpp'
hdr = hdr + ' Matrix.hpp MatrixImpl.hpp MatrixIO.hpp MatrixOrder.hpp MatrixRead.hpp'
hdr = hdr + ' MatrixStorage.hpp MatrixUtils.hpp MatrixWrite.hpp ParserDriver.hpp'
hdr = hdr + ' Parser.hpp pdb.hpp pdb_remarks.hpp pdbtraj.hpp PeriodicBox.hpp psf.hpp'
hdr = hdr + ' Selectors.hpp sfactories.hpp StreamWrapper.hpp timer.hpp'
hdr = hdr + ' TimeSeries.hpp tinker_arc.hpp tinkerxyz.hpp Trajectory.hpp'
hdr = hdr + ' UniqueStrings.hpp utils.hpp XForm.hpp'

loos_hdr_inst = env.Install(PREFIX + '/include', Split(hdr))

env.Alias('lib_install', [loos_lib_inst, loos_hdr_inst])

Return('loos')
