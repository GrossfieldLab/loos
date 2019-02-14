#!/usr/bin/env python3

import loos
import numpy

try:
    from scipy.spatial import ConvexHull as CH
except ImportError:
    import sys
    print("""
    Failure importing scipy, needed for ConvexHull calculations
    You either need to install scipy or adjust your PYTHONPATH
    to point to it.
    """)
    sys.exit(1)


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


def test_side(p, p0, p1):
    """
    Return val < 0 if p is to the left of the p0-p1 line segment,
           val > 0 if p is to the right
           val = 0 if p is on the line
    """
    val = (p.y() - p0.y()) * (p1.x() - p0.x()) - \
        (p.x() - p0.x()) * (p1.y() - p0.y())
    return val


class ConvexHull:

    def __init__(self, atomicgroup):
        self.atoms = atomicgroup
        self.hull = None
        self.vertices = None
        self.simplices = None

    def num_atoms(self):
        return len(self.atoms)

    def update_atoms(self, atomicgroup):
        self.atoms = atomicgroup

    def generate_hull(self):
        points = []
        for a in self.atoms:
            points.append([a.coords().x(), a.coords().y()])

        points = numpy.array(points)
        self.hull = CH(points)

    def generate_vertices(self):
        n = self.hull.simplices
        sorted_neighbors = []
        sorted_neighbors.append(n[0][0])
        sorted_neighbors.append(n[0][1])
        prev_entered = n[0][0]
        last_entered = n[0][1]
        num_simplices = len(self.hull.simplices)
        while (len(sorted_neighbors) < num_simplices):
            # search neighbor list to find the one with the current simplex
            # add its pair to the list
            last_entered = sorted_neighbors[-1]
            prev_entered = sorted_neighbors[-2]
            for i in range(len(n)):
                if ((n[i][0] == last_entered) and (n[i][1] != prev_entered)):
                    sorted_neighbors.append(n[i][1])
                    break
                if ((n[i][1] == last_entered) and (n[i][0] != prev_entered)):
                    sorted_neighbors.append(n[i][0])
                    break
        self.vertices = sorted_neighbors

    def atom(self, index):
        return self.atoms[int(index)]

    def coords(self, index):
        return self.atoms[int(index)].coords()

    def is_inside(self, p):
        """
        Returns true if p (a GCoord) is inside the hull
        """
        match = True
        side = test_side(p, self.coords(self.vertices[0]),
                         self.coords(self.vertices[1]))
        prev_side = side
        for i in range(1, len(self.vertices) - 1):
            side = test_side(p, self.coords(self.vertices[i]),
                             self.coords(self.vertices[i + 1]))
            match = (((side >= 0) and (prev_side >= 0)) or
                     ((side <= 0) and (prev_side <= 0)))
            if not match:
                return False
            prev_side = side

        # test the closing line segment
        side = test_side(p, self.coords(self.vertices[-1]),
                         self.coords(self.vertices[0]))
        match = (((side >= 0) and (prev_side >= 0)) or
                 ((side <= 0) and (prev_side <= 0)))
        return match


if __name__ == '__main__':

    ag = loos.createSystem("rhod_only.pdb")
    slicer = ZSliceSelector(-3.0, 3.0)

    ag = loos.selectAtoms(ag, 'name == "CA"')
    ag_slice = slicer(ag)

    hull = ConvexHull(ag_slice)
    hull.generate_hull()
    hull.generate_vertices()

    i = 0
    for v in hull.vertices:
        c = hull.coords(v)
        print(i, v, c.x(), c.y(), c.z())
        i += 1

    print(hull.is_inside(loos.GCoord(0.0, 0.0, 0.0)))
    print(hull.is_inside(loos.GCoord(20.0, 0.0, 0.0)))
    print(hull.is_inside(loos.GCoord(0.0, 20.0, 0.0)))
