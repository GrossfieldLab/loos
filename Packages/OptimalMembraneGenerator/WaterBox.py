#!/usr/bin/env python3

import loos

# @cond TOOLS_INTERNAL


class WaterBox:
    def __init__(self, filename, template_box, target_box, segname, num_sites):
        self.filename = filename
        self.box = target_box
        self.template_box = template_box
        self.segname = segname
        self.num_sites = num_sites

        self.template = loos.createSystem(filename)
        for i in range(len(self.template)):
            self.template[i].segid(segname)

        self.buildBox()

    def buildBox(self):
        """
        Use the template structure to build the full structure by
        replicating it each dimension until it's bigger than the
        target box size, then deleting waters with centroid that fall
        outside the target box.
        The resulting AtomicGroup is self.full_system
        """

        # figure out how many replicas we need in each direction
        num_x = int(self.box.x() / self.template_box.x()) + 1
        num_y = int(self.box.y() / self.template_box.y()) + 1
        num_z = int(self.box.z() / self.template_box.z()) + 1

        # build up the new system by replicating the target box
        # and systematically translating it
        self.full_system = loos.AtomicGroup()
        for x in range(num_x):
            for y in range(num_y):
                for z in range(num_z):
                    new = self.template.copy()

                    trans = loos.GCoord()
                    trans.x(self.template_box.x() * x)
                    trans.y(self.template_box.y() * y)
                    trans.z(self.template_box.z() * z)

                    new.translate(trans)
                    self.full_system.append(new)

        self.full_system.centerAtOrigin()

        # trim the waters outside the target box size
        residues = self.full_system.splitByResidue()
        half_box = 0.5 * self.box
        to_remove = loos.AtomicGroup()
        for res in residues:
            centroid = res.centroid()
            if ((abs(centroid.x()) > half_box.x()) or
                (abs(centroid.y()) > half_box.y()) or
                (abs(centroid.z()) > half_box.z())):
                    to_remove.append(res)

        print("Need to remove: ", len(to_remove))
        print("before: ", self.full_system.boundingBox(), len(self.full_system))
        self.full_system.remove(to_remove)
        print("after: ", self.full_system.boundingBox(), len(self.full_system))

        self.full_system.periodicBox(self.box)

        # renumber atom ids and resids
        self.full_system.renumber()
        for i in range(len(self.full_system)):
            self.full_system[i].resid(i//3 + 1)
        #residues = self.full_system.splitByResidue()
        #for i in range(len(residues)):
        #    for j in range(len(residues[i])):
        #        residues[i][j].resid(i+1)

    def append_waters(self, other):
        """
        "other" is an AtomicGroup of waters.  Merge them
        into the current WaterBox, renumbering the atoms
        and residues, and updating their segment name
        """

        self.full_system += other
        self.full_system.renumber()
        residues = self.full_system.splitByResidue()
        for i in range(len(residues)):
            for j in range(len(residues[i])):
                residues[i][j].resid(i+1)
                residues[i][j].segid(self.segname)

    def pdb(self):
        """
        Return a string containing a PDB version of the full_system,
        convenient for writing out the coordinates in PDB format.
        """

        p = loos.PDB.fromAtomicGroup(self.full_system)
        return str(p)

# @endcond


if __name__ == '__main__':
    coordfile = 'water_small.crd'
    box_size = loos.GCoord(15.5516, 15.5516, 15.5516)
    big_box = loos.GCoord(74.1, 74.1, 95.0)

    w = WaterBox(coordfile, box_size, big_box, "BULK")

    f = open("big_water.pdb", "w")
    f.write(w.pdb())
    f.close()

    print(w.full_system.periodicBox())
    print(w.full_system.boundingBox())
