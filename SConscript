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
apps = apps + ' ccpdb.cpp pdbtraj.cpp tinker_arc.cpp KernelStack.cpp'


loos = env.Library('loos', Split(apps))

Return('loos')
