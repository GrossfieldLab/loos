#!/usr/bin/env python3
"""
NAMDBin reads and writes the binary format used by NAMDBin
for restart files.  I've written it to support coordinates
and velocities, because that's what LOOS supports, but the
same format is also used to store forces.

Alan Grossfield
University of Rochester Medical Center
2017
"""
import struct
import loos

class NAMDBin:
    def __init__(self, filename, model=None):
        self.filename = filename
        self.atomicGroup = model
        self.endian_signal = ">"  # supply a default value
        # NAMD stores velocities in weird unit, this converts
        # to Ang/ps
        self.PDBVELFACTOR = 20.45482706

    def set_little_endian(self):
        self.endian_signal = "<"  # supply a default value

    def set_big_endian(self):
        self.endian_signal = ">"  # supply a default value

    def readVals(self, velocities=False):
        """
        Read the binary format.  Determine endian-ness by
        verifying the number of atoms (the first thing read)
        matches the supplied AtomicGroup.
        """
        if self.atomicGroup is None:
            raise ValueError("Must supply an atomic group first")

        with open(self.filename, 'rb') as f:
            atoms_bytes = f.read(4)
            # check the number of atoms to figure out endian-ness
            num_atoms_big = struct.unpack('>i', atoms_bytes)[0]
            num_atoms_little = struct.unpack('<i', atoms_bytes)[0]
            num_atoms_model = len(self.atomicGroup)
            if num_atoms_big == num_atoms_model:
                self.endian_signal = ">"
            elif num_atoms_little == num_atoms_model:
                self.endian_signal = "<"
            else:
                raise ValueError("Number of atoms doesn't match supplied model")

            for i in range(num_atoms_model):
                x, y, z = struct.unpack(self.endian_signal + 'ddd', f.read(24))
                c = loos.GCoord(x, y, z)
                if velocities:
                    c *= self.PDBVELFACTOR
                    self.atomicGroup[i].velocities(c)
                else:
                    self.atomicGroup[i].coords(c)

    def writeVals(self, filename, velocities=False):
        """
        Write a binary coordinate/velocity file, using the
        same endian-ness as our input.
        """
        with open(filename, 'wb') as f:
            bs = struct.pack(self.endian_signal + 'i', len(self.atomicGroup))
            f.write(bs)
            for i in range(len(self.atomicGroup)):
                if velocities:
                    c = self.atomicGroup[i].velocities() / self.PDBVELFACTOR
                else:
                    c = self.atomicGroup[i].coords()
                bs = struct.pack(self.endian_signal + 'ddd',
                                 c.x(), c.y(), c.z())
                f.write(bs)


if __name__ == '__main__':
    # Testing pieces
    from copy import copy

    model_file = "namd.psf"
    model = loos.createSystem(model_file)
    model2 = copy(model)
    nb = NAMDBin("npgt_restart.100.vel", model)
    nb.readVals(velocities=True)
    nb.writeVals("foo.vel", velocities=True)

    nb2 = NAMDBin("foo.vel", model2)
    nb2.readVals(velocities=True)
    for i in range(len(model)):
        m1 = model[i].velocities()
        m2 = model2[i].velocities()
        d = (m2 - m1).length()

        print(model[i].velocities(), model2[i].velocities(), d)
