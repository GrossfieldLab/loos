from loos import *

mol = createSystem("lfb.psf")
traj = createTrajectory("lfb.dcd", mol)
backbone = selectAtoms(mol, "name =~ '^(C|O|N|CA)$' && segid == 'PE1'")

t = 0
while (traj.readFrame()):
    traj.updateGroupCoords(mol)
    axes = backbone.principalAxes()
    print t, "\t", abs(axes[0].z())
    t += 1


