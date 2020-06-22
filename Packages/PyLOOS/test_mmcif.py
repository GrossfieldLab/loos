#!/usr/bin/env python3

import loos

#mmcif_filename = "/Users/agrossfield/Downloads/1u19.cif"
#mmcif = loos.MMCIF(mmcif_filename)

#print(mmcif.unitCell())
#print(mmcif.periodicBox())



#mmcif2 = loos.createSystem(mmcif_filename)
#print(mmcif2[0])
#
#
#sys3 = loos.createSystem("/Users/agrossfield/Downloads/1u19.pdbx")
#print(type(sys3))

sys3 = loos.createSystem("/Users/agrossfield/Downloads/1u19.pdb")

mmcif2 = loos.MMCIF.fromAtomicGroup(sys3)
print(type(mmcif2))

#print(str(mmcif))

with open("/Users/agrossfield/Downloads/new.cif", "w") as outfile:
    outfile.write(str(mmcif2))

print(str(mmcif2))
