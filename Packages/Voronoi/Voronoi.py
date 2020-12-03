#!/usr/bin/env python3
#  Do voronoi analysis on slices of a membrane
#  2D only
#  This code is a wrapper around the scipy Voronoi class, which is itself
#  a wrapper around QHull (http://www.qhull.org/)
#
#  Alan Grossfield, University of Rochester Medical Center, Dept
#                   of Biochemistry and Biophysics
#  Copyright 2014
#

import sys
import loos
try:
    import numpy
except ImportError as e:
    print("Error importing numpy: ({0}): {1}".format(e.errno, e.strerror))
try:
    from scipy.spatial import Voronoi
except ImportError as e:
    print("Error importing Voronoi from scipy: ({0}): {1}".format(
        e.errno, e.strerror))


class ZSliceSelector:
    def __init__(self, min_z, max_z):
        self.min_z = float(min_z)
        self.max_z = float(max_z)

    def __call__(self, atomicgroup):
        ag = loos.AtomicGroup()
        ag.periodicBox(atomicgroup.periodicBox())
        for a in atomicgroup:
            if self.min_z < a.coords().z() < self.max_z:
                ag.append(a)
        return ag


class VoronoiError(Exception):
    """Base class for Voronoi package exceptions"""

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return(repr(self.msg))


class VoronoiWrapper:
    """
    Wrap the scipy Voronoi class, which in turn is a wrapper around QHull.
    So, it's a wrapper wrapper, and if this program is called from a script,
    then we'd have a wrapper wrapper wrapper.  Now I just need m4 and Tod will
    approve.

    atomicgroup is a LOOS AtomicGroup
    pad is a float specifying how far out we will generate padding atoms
        (fake "image" atoms used to emulate periodicity).  15 ang is a good
        default if you're using all atoms or all heavy atoms, but you may
        need to go farther if you're using a sparser selection (e.g. just
        lipid phosphates)
    """

    def __init__(self, atomicGroup, pad=15.0):
        self.atoms = atomicGroup
        self.edges = []
        self.pad = pad
        self.padding_atoms = []
        self.voronoi = None
        self.regions = []
        self.superRegions = []
        self.atoms_to_regions = {}

    def isPeriodic(self):
        return self.atoms.isPeriodic()

    def num_atoms(self):
        return len(self.atoms)

    def update_atoms(self, atomicGroup):
        self.atoms = atomicGroup

    def generate_padding_atoms(self):
        """
        Build the list of atoms in the periodic images surrounding the central
        image.

        This is necessary because QHull doesn't know anything about periodicity,
        so we need to do fake it; otherwise, you'd have bizarre or infinite
        areas for atoms near the edge of the box
        """
        if not self.isPeriodic():
            raise VoronoiError("Periodic boundaries are required")

        box = self.atoms.periodicBox()
        half_box = 0.5 * box
        xbox = loos.GCoord(box.x(), 0.0, 0.0)
        ybox = loos.GCoord(0.0, box.y(), 0.0)

        for a in self.atoms:
            c = a.coords()
            # generate the square neighbors
            if (c.x() < -half_box.x() + self.pad):
                self.padding_atoms.append(c + xbox)
            if (c.x() > half_box.x() - self.pad):
                self.padding_atoms.append(c - xbox)
            if (c.y() < -half_box.y() + self.pad):
                self.padding_atoms.append(c + ybox)
            if (c.y() > half_box.y() - self.pad):
                self.padding_atoms.append(c - ybox)

            # generate the diagonal neighbors
            if ((c.x() < -half_box.x() + self.pad) and
                    (c.y() < -half_box.y() + self.pad)):
                self.padding_atoms.append(c + xbox + ybox)

            if ((c.x() > half_box.x() - self.pad) and
                    (c.y() > half_box.y() - self.pad)):
                self.padding_atoms.append(c - xbox - ybox)

            if ((c.x() < -half_box.x() + self.pad) and
                    (c.y() > half_box.y() - self.pad)):
                self.padding_atoms.append(c + xbox - ybox)

            if ((c.x() > half_box.x() - self.pad) and
                    (c.y() < -half_box.y() + self.pad)):
                self.padding_atoms.append(c - xbox + ybox)

        return len(self.padding_atoms)

    def num_padding_atoms(self):
        return len(self.padding_atoms)

    def get_region_from_atomid(self, atomid):
        return self.atoms_to_regions[atomid]

    def generate_voronoi(self):
        # make the 2D list of points
        points = []
        for a in self.atoms:
            points.append([a.coords().x(), a.coords().y()])

        numpad = self.generate_padding_atoms()
        for p in self.padding_atoms:
            points.append([p.x(), p.y()])

        points = numpy.array(points)
        self.voronoi = Voronoi(points)

        # read in all of the ridges
        for r in self.voronoi.ridge_vertices:
            self.edges.append(Edge(r[0], r[1]))

        # now build the regions
        # the scipy voronoi object has an array point_region, which
        # maps the index of the input point to its associated Voronoi region
        v = []
        for vert in self.voronoi.vertices:
            c = loos.GCoord(vert[0], vert[1], 0.0)
            v.append(c)
        for i in range(self.num_atoms()):
            index = self.voronoi.point_region[i]
            if index == -1:
                raise ValueError("point %d (atomId = %d) from voronoi decomposition isn't associated with a voronoi region; you may need to increase the padding value" % i).with_traceback(
                    self.atoms[i].id())
            r = self.voronoi.regions[index]
            self.regions.append(Region(v, r, self.atoms[i]))
            self.atoms_to_regions[self.atoms[i].id()] = self.regions[i]


class Edge:
    def __init__(self, index1, index2):
        self.ind1 = index1
        self.ind2 = index2

    def __eq__(self, other):
        if ((self.ind1 == other.ind1) and (self.ind2 == other.ind2) or
                (self.ind2 == other.ind1) and (self.ind1 == other.ind2)):
            return True
        else:
            return False


class Region:
    def __init__(self, vert_array, indices, atom):
        if len(indices) == 0:
            raise ValueError("can't have 0-length list of indices")
        self.vertices = vert_array
        self.indices = indices
        self.atom = atom
        self.edges = []

        # build list of indices
        for i in range(len(indices) - 1):
            self.edges.append(Edge(indices[i], indices[i + 1]))
        self.edges.append(Edge(indices[-1], indices[0]))

    def num_indices(self):
        return len(self.indices)

    def num_edges(self):
        return len(self.edges)

    def area(self):
        area = 0.0
        for i in range(self.num_indices()):
            if ((i + 1) == self.num_indices()):
                j = 0
            else:
                j = i + 1

            p1 = self.vertices[self.indices[i]]
            p2 = self.vertices[self.indices[j]]
            area += p1.x() * p2.y() - p2.x() * p1.y()

        area = 0.5 * abs(area)
        return(area)

    def is_neighbor(self, other):
        """
        Two regions are neighbors if they have an edge in common
        """
        for e1 in self.edges:
            for e2 in other.edges:
                if (e1 == e2):
                    return(True)
        return(False)

    def print_indices(self):
        for i in range(self.num_indices()):
            print(self.indices[i], " ", self.vertices[self.indices[i]])

    def atomId(self):
        return self.atom.id()

    def coords(self):
        return self.atom.coords()

    def __str__(self):
        string = ""
        for i in range(self.num_indices()):
            g = self.vertices[self.indices[i]]
            string += "%f\t%f\n" % (g.x(), g.y())
        g = self.vertices[self.indices[0]]
        string += "%f\t%f\n" % (g.x(), g.y())
        return string


class SuperRegion:

    def __init__(self, regions=None):
        if regions:
            self.regions = regions
        else:
            self.regions = []

    def add_region(self, region):
        self.regions.append(region)

    def buildFromAtoms(self, atomicGroup, voronoi):
        for atom in atomicGroup:
            self.add_region(voronoi.get_region_from_atomid(atom.id()))

    def area(self):
        a = 0.0
        for r in self.regions:
            a += r.area()
        return(a)

    def is_neighborRegion(self, other):
        for r in self.regions:
            if r.is_neighbor(other):
                return True
        return False

    def is_neighborSuperRegion(self, other):
        for r in self.regions:
            for r2 in other.regions:
                if r.is_neighbor(r2):
                    return True
        return False

    def print_indices(self):
        for r in self.regions:
            r.print_indices()
            print()


if __name__ == '__main__':

    import sys

    structure = loos.createSystem("trj_1.pdb")
    #structure = loos.createSystem("example.pdb")
    #structure = loos.createSystem("b2ar.pdb")

    #box = loos.GCoord(55., 77, 100)
    #box = loos.GCoord(55., 77, 100)
    phos = loos.selectAtoms(structure, 'name == "P"')
    upper = loos.AtomicGroup()
    for p in phos:
        if p.coords().z() > 0:
            upper.append(p)
    upper.periodicBox(structure.periodicBox())
    print(upper.isPeriodic())
    """
    slice = loos.AtomicGroup()
    for a in structure:
        if a.coords().z() > 20 and a.coords().z() < 21:
            slice.append(a)
    slice.periodicBox(box)
    """

    v = VoronoiWrapper(upper)
    #v = VoronoiWrapper(slice)
    # print v.isPeriodic()
    # print v.num_atoms()
    v.generate_voronoi()
    print(v.num_atoms(), v.num_padding_atoms())
    for i in range(v.num_atoms()):
        print(i, v.atoms[i].coords())
        print(v.regions[i].area())

    for i in range(v.num_padding_atoms()):
        print(i, v.padding_atoms[i])

    s = SuperRegion(v.regions[:1])
    s2 = SuperRegion(v.regions[:2])
    print(s.area(), v.regions[0].area())
    print(s2.area(), v.regions[0].area() + v.regions[1].area())
    s.add_region(v.regions[1])
    print(s.area())
    s.print_indices()

    s3 = SuperRegion()
    s3.buildFromAtoms(upper, v)
    print(s3.area())
    total = 0.0
    for r in v.regions:
        total += r.area()
        print(r.coords())
        print(r)
    print(total)
