#!/usr/bin/env python3

#
#  This file is part of LOOS.
#
#  LOOS (Lightweight Object-Oriented Structure library)
#  Copyright (c) 2012 Tod D. Romo, Grossfield Lab
#  Department of Biochemistry and Biophysics
#  School of Medicine & Dentistry, University of Rochester
#
#  This package (LOOS) is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation under version 3 of the License.
#
#  This package is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy


def read_file(filename):
    file = open(filename)
    v = []
    z_vals = []
    comments = []
    for line in file.xreadlines():
        if line.startswith("#"):
            comments.append(line)
            continue
        fields = line.split()
        z = float(fields[0])
        temp = map(float, fields[1:])
        v.append(temp)
        z_vals.append(z)
    v = numpy.array(v)
    z_vals = numpy.array(z_vals)
    return z_vals, v, comments


if __name__ == '__main__':

    import sys

    if (len(sys.argv) > 1 and sys.argv[1] == "--fullhelp"):
        print("""
SYNOPSIS

Compute electrostatic potential along membrane normal

DESCRIPTION

This program generates the electrostatic potential profile for a membrane
system given a data file containing the charge density as a function of position
along the membrane normal.  Intended to post-process the output of the
density-dist tool, it takes input in electrons/Ang^3 and outputs the potential
in volts.

The algorithm used here is described by Sachs, et al, J Chem Phys, 2004, 121,
10847.  In particular, the periodicity correction is applied, such that the top
and bottom of the periodic box must have the same electrostatic potential.
However, this correction is applied _only_ to the full charge density, not the
charge density of each individual component (assuming the user made selections
when running density-dist).  This is because only the total potential needs to
be continuous at the potential, not that due to any given component.

EXAMPLE

potential_profile.py is intended to be used in combination with the LOOS
tool density-dist.  For example:

density-dist --type=charge -- path/to/model-file path/to/trajectory -38 38 76 > charge-density.dat
potential_profile.py charge-density.dat > potential.dat

The "--" in the first line tells LOOS to stop interpreting arguments with
leading minus signs as flags; otherwise, it will choke on the "-38" (run
density-dist --fullhelp for more discussion).

The first column of the output is the z-value, the second is the periodicity
correction (see above), and the subsequent columns are the electrostatic
potentials for the full system and any individual components selected when
density-dist was run.

              """)
        sys.exit(1)
    elif (len(sys.argv) != 2):
        print("Usage: ", sys.argv[0], " filename")
        print("  where filename is presumed to be output from density-dist")
        sys.exit(1)

    datafilename = sys.argv[1]
    z_vals, data, comments = read_file(datafilename)

    # assume input units are angstroms and electron charge,
    # convert to V
    eps0 = 8.85e-12    # C/V m
    e_to_C = 1.60e-19  # C/e
    ang_to_m = 1e-10

    #      dielectric   charge density        integrated distance
    units = 1/eps0 * e_to_C/(ang_to_m**3) * ang_to_m**2

    pot = numpy.zeros(data.shape, float)
    # compute the potential, uncorrected for periodicity
    for i in range(len(data)):
        for j in range(i+1, len(data)):
            dz = z_vals[i] - z_vals[j]
            pot[j] += dz*data[i]

    # compute the periodicity correction for the first column only!!!
    zrange = z_vals[-1] - z_vals[0]
    corr = z_vals.copy()
    delta_p = pot[-1][0] - pot[0][0]
    corr *= -delta_p/zrange

    pot = pot.swapaxes(0, 1)
    pot[0] += corr
    # try recentering in the middle: can help if there are artifacts at the
    # upper edge
    pot[0] -= pot[0][len(pot[0])//2]

    pot = pot.swapaxes(0, 1)

    # convert to units of volts
    pot *= units
    corr *= units

    # output the result
    print("# ", " ".join(sys.argv))
    print("#Generated from: ")
    print("#".join(comments),)
    print("#z\tCorrection\tPotentials")

    for i in range(len(z_vals)):
        s = " ".join(map(str, pot[i].tolist()))
        print(z_vals[i], corr[i], s)
