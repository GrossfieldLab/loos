#!/usr/bin/env python3

import subprocess
import sys
import os.path
import loos


# @cond TOOLS_INTERNAL
class NAMD:

    def __init__(self, psf_file, start_pdb, end_pdb, param_file, box,
                 command):
        self.psf_file = psf_file
        self.start_pdb = start_pdb
        self.end_pdb = end_pdb
        self.param_file = param_file

        if command is None:
            self.command = "/opt/bin/namd2"
        else:
            self.command = command

        # in case we're going to do constraints, construct this filename
        self.cons_k_filename = self.start_pdb[:-4] + ".cons.pdb"

        self.x_box = box.x()
        self.y_box = box.y()
        self.z_box = box.z()

    def update_box(self, box):
        self.x_box = box.x()
        self.y_box = box.y()
        self.z_box = box.z()

    def construct_header(self):
        lines = []
        lines.append("structure " + self.psf_file)
        lines.append("coordinates " + self.start_pdb)
        lines.append("outputname " + self.end_pdb)
        lines.append("paratypecharmm on")
        for p in self.param_file:
            lines.append("parameters " + p)
        lines.append("exclude scaled1-4")
        lines.append("1-4scaling 1.0 ")
        lines.append("binaryoutput no")

        string = "\n".join(lines)
        string += self.template() + "\n"
        return string

    def construct_box(self):
        lines = []
        lines.append("cellBasisVector1 " + str(self.x_box) + " 0 0")
        lines.append("cellBasisVector2 0 " + str(self.y_box) + " 0")
        lines.append("cellBasisVector3 0 0 " + str(self.z_box))

        string = "\n".join(lines)
        string += "\n"
        return string

    def construct_mini(self, num_iter=100):
        line = "minimize " + str(num_iter) + "\n"
        return line

    def construct_dyn(self, num_iter=100):
        line = "run " + str(num_iter) + "\n"
        return line

    def write_inputfile(self, filename, nsteps=100):
        file = open(filename, "w")
        file.write(self.construct_header())
        file.write(self.construct_constraints())
        file.write(self.construct_box())
        file.write(self.construct_mini(nsteps))
        file.write(self.construct_dyn(nsteps))
        file.close()

    def write_restraintfile(self, directory, atomicgroup, spring=10.0):
        pdb = loos.PDB.fromAtomicGroup(atomicgroup.copy())
        spring_coord = loos.GCoord(spring, 0., 0.)
        for atom in pdb:
            atom.coords(loos.GCoord(0., 0., 0.))
        heavy = loos.selectAtoms(pdb, '!(name =~ "^H")')
        for atom in heavy:
            atom.coords(spring_coord)

        pdb_file = open(os.path.join(directory, self.cons_k_filename), "w")
        pdb_file.write(str(pdb))
        pdb_file.close()

    def construct_constraints(self):
        lines = [ "constraints on",
                  "selectConstraints on",
                  "selectConstrZ on",
                  "conskcol X"
                ]
        line = "consref " + self.start_pdb
        lines.append(line)

        line = "conskfile " + self.cons_k_filename
        lines.append(line)

        line = "\n".join(lines)
        line += "\n"
        return line

    def run_namd(self, inputfilename, outfilename):
        """
        Run namd on 4 processors, and report any failure to stderr
        """
        outfile = open(outfilename, "w")
        try:
            subprocess.check_call([self.command, "+p4", inputfilename],
                                  stdout=outfile)
        except subprocess.CalledProcessError:
            sys.stderr.write("NAMD call failed, inp = %s, out = %s\n" %
                             (inputfilename, outfilename))
            sys.exit(-1)

    def template(self):
        """
        Template entries for the NAMD input file.  Note that these are
        note settings you'd want to use for any kind of production run
        (particularly the insanely short cutoff).  However, we're mainly
        using these runs for local minimization and clash-alleviation, so
        this approach is fine and much faster than doing a realistic cutoff
        and Ewald summation.
        """
        return """
temperature 500
stepsPerCycle  20
switching on
switchdist 5
cutoff 6
pairlistdist 7
COMmotion no
zeroMomentum yes
rigidBonds all
rigidTolerance 0.0000000001
useSettle on
wrapAll on
outputenergies 50
timestep 2.0
langevin on
langevinTemp 500
langevinDamping 2
langevinHydrogen off

"""
# @endcond


if __name__ == '__main__':
    import loos

    box = loos.GCoord(377, 377, 1000)
    n = NAMD("generated.psf", "t.pdb", "end", "toppar/par_build.inp", box)
    print(n.construct_header())

    print(n.construct_box())

    box[0] = 60
    n.update_box(box)
    print(n.construct_box())
    box[0] = 377
    n.update_box(box)

    print(n.construct_mini())

    n.write_inputfile("n.inp", 50)
    print("launching NAMD")
    n.run_namd("n.inp", "n.out")
    print("finished NAMD")
