#!/usr/bin/env python3
import sys
import loos
import loos.OptimalMembraneGenerator
from loos.OptimalMembraneGenerator import LipidLibrary
import subprocess
import os

# @cond TOOLS_INTERNAL


class Segment:
    """
    Class describing a segment of a PSF file.  Contains the information
    needed to write a psf, and to build this set ouf lipids.

    The input value line is expected to look like:
SEGMENT TPC       POPC     90    19      P   1      ./popc_c36
    where:
        TPC =  segment name (all of 1 lipid type, 1 res/lipid)
        POPC=  residue name
        90  =  number of lipids (and residues) in the segment
        19  =  distance of phosphate from center of bilayer
        P   =  atom to be placed at the distance (usually P for lipid phosphate)
        1   =  placement (1 for upper leaflet, -1 for lower, 0 for TM)
  ./popc_c36=  path to the lipid library where the configurations are stored
    """
    def __init__(self, line):
        (segment, segname, resname, number, phos_pos, phos_atom, placement, library) = line.split()

        self.segname = segname
        self.resname = resname
        self.numres = int(number)
        self.library = LipidLibrary.LipidLibrary(library)

        self.placement = int(placement)
        self.phos_height = float(phos_pos)
        self.phos_atom = phos_atom

    def write(self):
        arr = ["segment " + self.segname + " {"]
        for i in range(1, self.numres + 1):
            arr.append("    residue " + str(i) + " " + self.resname)

        arr.append("}")
        string = "\n".join(arr)

        return string


class WaterSeg(Segment):
    """
    Subclass of Segment (done to inherit the write() method) intended to handle
    water.
    """
    def __init__(self, line):
        fields = line.split()
        self.segname = fields[1]
        self.resname = fields[2]
        self.numres = int(fields[3])
        self.thickness = float(fields[4])
        self.box_size = float(fields[5])
        self.coords_filename = fields[6]
        if len(fields) > 7:
            self.num_sites = int(fields[7])
        else:
            self.num_sites = 3


    def write(self):
        arr = ["segment " + self.segname + " {"]
        arr.append("    auto none")
        for i in range(1, self.numres + 1):
            arr.append("    residue " + str(i) + " " + self.resname)

        arr.append("}")
        string = "\n".join(arr)

        return string


class SaltSeg(Segment):
    """
    Subclass of Segment (done to inherit the write() method) intended to handle
    single atom ions.
    """
    def __init__(self, line):
        (tag, segname, resname, number) = line.split()
        self.segname = segname
        self.resname = resname
        self.atomname = resname  # assume atom name is the same as the residue
        self.numres = int(number)


class Protein:
    """
    Analogous to Segment class, handles things like proteins embedded in the
    membrane
    """
    def __init__(self, line):
        (tag, model_file, psf_file, water_seg, scale_by_molecule) = line.split()
        self.model_file = os.path.abspath(model_file)
        self.psf_file = os.path.abspath(psf_file)
        self.water_segname = water_seg.upper()
        self.scale = int(scale_by_molecule)

        self.model = loos.createSystem(self.model_file)
        self.segments = self.model.splitByUniqueSegid()

        self.has_water = False
        for s in self.segments:
            if s[0].segid() == self.water_segname:
                self.has_water = True
                break

    def is_water(self, segname):
        return segname == self.water_segname

    def water_seg(self):
        for s in self.segments:
            if s[0].segid() == self.water_segname:
                return s
        return None

    def write(self, include_water=False):
        """
        Write psfgen input for all of the segments EXCEPT the water
        segment
        """

        string = "readpsf " + self.psf_file + "\n"
        if not include_water:
            string += "# remove the water segment from the psf because\n"
            string += "# waters are part of the main water segment now\n"
            string += "delatom " + self.water_segname + "\n\n"

        return string


class ReadConfig:
    """
    Class to read the config file and set up the segments.  Used to drive
    psfgen and to guide lipid construction.
    """
    def __init__(self, filename):

        self.segments = []
        self.water = None
        self.salt = []
        self.protein = None
        self.protrot = None

        self.topology = []
        self.parameters = []
        self.psfname = None
        self.directory = "./out/"

        self.box = None

        self.namd_binary = None
        self.psfgen_binary = None

        file = open(filename)

        for line in file.readlines():
            # skip blanks and comments
            if line.startswith("#") or line.isspace() or len(line) == 0:
                continue
            if line.upper().startswith("TOPOLOGY"):
                (top, filename) = line.split()
                self.topology.append(os.path.abspath(filename))
            elif line.upper().startswith("PARAMETERS"):
                (par, filename) = line.split()
                self.parameters.append(os.path.abspath(filename))
            elif line.upper().startswith("PSFGEN"):
                (n, psfgen) = line.split()
                self.psfgen_binary = psfgen
            elif line.upper().startswith("PSF"):
                (p, psfname) = line.split()
                self.psfname = psfname
            elif line.upper().startswith("SEGMENT"):
                s = Segment(line)
                self.segments.append(s)
            elif line.upper().startswith("WATER"):
                self.water = WaterSeg(line)
            elif line.upper().startswith("SALT"):
                s = SaltSeg(line)
                self.salt.append(s)
            elif line.upper().startswith("PROTEIN"):
                self.protein = Protein(line)
            elif line.upper().startswith("PROTROT"):
                (r, protrot) = line.split()
                self.protrot = int(protrot)
            elif line.upper().startswith("BOX"):
                (b, x, y, z) = line.split()
                x = float(x)
                y = float(y)
                z = float(z)
                self.box = loos.GCoord(x, y, z)
            elif line.upper().startswith("NAMD"):
                (n, namd) = line.split()
                self.namd_binary = namd

            else:
                sys.stderr.write("Unrecognized line type: %s" % line)

        if len(self.topology) == 0:
            sys.stderr.write("No topology file specified... exiting\n")
            sys.exit(1)

        if len(self.parameters) == 0:
            sys.stderr.write("No parameter file specified... exiting\n")
            sys.exit(1)

        if self.psfname is None:
            sys.stderr.write("No output psf file specified... exiting\n")
            sys.exit(1)

        # Warn but don't exit if there are no lipids
        if len(self.segments) == 0:
            sys.stderr.write("No segments specified... your system has no lipids\n")

        # TODO: may want to add similar checks for the other keys

    def generate_psf(self,
                     include_lipid=True,
                     include_water=False,
                     include_other=False,
                     include_other_water=True,
                     output_psf=None):

        if output_psf is None:
            output_psf = self.psfname

        lines = []
        for t in self.topology:
            line = "topology " + t
            lines.append(line)

        lines.append("\n")

        if include_other and self.protein is not None:
            line = self.protein.write(include_other_water)
            lines.append(line)

        lines.append("\n")

        if include_lipid:
            for s in self.segments:
                line = s.write()
                lines.append(line)
                lines.append("\n")

        line = ""
        if include_other:
            if self.water is not None and include_water:
                line = self.water.write()
                lines.append(line)
                if len(self.salt) > 0:
                    for s in self.salt:
                        line = s.write()
                        lines.append(line)
                lines.append("\n")

        lines.append("writepsf x-plor cmap " + output_psf)
        a = "\n".join(lines)
        return a

    def __getitem__(self, index):
        return self.segments[index]


class PSFGen:
    def __init__(self, psf_string, command=None):
        if command is None:
            self.command = "/opt/bin/psfgen"
        else:
            self.command = command

        self.psf_string = psf_string

    def run(self):
        psfgen = subprocess.Popen("", 0,
                                  self.command,
                                  stdin=subprocess.PIPE,
                                  universal_newlines=True)
        psfgen.communicate(self.psf_string)

# @endcond


if __name__ == '__main__':

    config_file = sys.argv[1]

    conf = ReadConfig(config_file)

    #print conf.generate_psf()
    print(conf.generate_psf(True, True))
    #print conf.generate_psf(True, False)
    #print conf.generate_psf(False, True)
