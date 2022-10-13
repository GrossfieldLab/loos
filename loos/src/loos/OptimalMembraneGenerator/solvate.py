#!/usr/bin/env python3

import sys
import os
import random
import stat
import shutil
import loos
import loos.OptimalMembraneGenerator
from loos.OptimalMembraneGenerator import PSFGen
from loos.OptimalMembraneGenerator import WaterBox

def main():
    config = PSFGen.ReadConfig(sys.argv[1])
    if len(sys.argv) > 2:
        config.directory = sys.argv[2]


    # make sure the working directory exists and is writeable
    if not os.path.exists(config.directory):
        sys.stderr.write("Creating directory %s\n" % config.directory)
        os.makedirs(config.directory)
    elif not os.access(config.directory, os.R_OK | os.W_OK | os.X_OK):
        sys.stderr.write("Don't have access to directory %s\n" % config.directory)
        sys.stderr.write("Will attempt to change permissions\n")
        try:
            os.chmod(config.directory, stat.S_IREAD | stat.S_IWRITE | stat.S_IEXEC)
        except OSError:
            sys.stderr.write(" Unable to change permissions on directory\n")
            os.exit(-1)


    # We need a model containing just the protein and lipid, so we'll
    # make a psfgen script, run it, use the psf to make the AtomicGroup,
    temporary_psfname = os.path.join(config.directory, config.psfname)
    psfgen_script = config.generate_psf(True, False, True, True, temporary_psfname)
    psfgen = PSFGen.PSFGen(psfgen_script, config.psfgen_binary)
    psfgen.run()
    system = loos.createSystem(temporary_psfname)

    segments = []
    for segment in config.segments:
        s = loos.selectAtoms(system, 'segname == "' + segment.segname + '"')
        if len(s) == 0:
            sys.stderr.write(
                "Selection failed assembling system: segment %s doesn't exist\n"
                % (segment.segname)
            )
            sys.stderr.write("Exiting...\n")
            sys.exit(0)

        segments.append(s)


    # copy the protein coordinates into the system
    if config.protein is not None:
        # create AtomicGroup containing all protein segments in case
        # we want to rotate it
        to_rot = loos.AtomicGroup()

        for s in config.protein.segments:
            current_seg = s[0].segid()
            # Don't need to trap failed selection here, because we
            # already know this segment exists
            seg = loos.selectAtoms(system, 'segname == "' + current_seg + '"')
            seg.copyMappedCoordinatesFrom(s)
            to_rot.append(seg)

        # if we asked for rotation, rotate all segments together
        # about a random axis
        if config.protrot:
            axis = loos.GCoord()
            axis.random()
            rot = random.uniform(0.0, 360.0)
            to_rot.rotate(axis, rot)


    sys.stderr.write("Beginning water box construction\n")
    # now add water and salt
    water_template = loos.GCoord(1.0, 1.0, 1.0)
    water_template *= config.water.box_size
    water_target = loos.GCoord(config.box.x(), config.box.y(), config.box.z())
    water = WaterBox.WaterBox(config.water.coords_filename,
                            water_template,
                            water_target,
                            config.water.segname,
                            config.water.num_sites
                            )


    # Figure out how many ions we're going to add
    total_salt = 0
    for salt in config.salt:
        total_salt += salt.numres
    total_water_and_salt = total_salt + config.water.numres

    sys.stderr.write(
        "Water box has %d waters before superposition\n"
        % (len(water.full_system) // water.num_sites)
    )
    sys.stderr.write("Final target: %d waters\n" % (config.water.numres))


    # Verify we have enough water.  We need enough to end up with
    # the planned number of waters, even after we remove one water molecule
    # for each ion we add.
    if len(water.full_system) // water.num_sites < total_water_and_salt:
        raise ValueError(
            "Too few waters before superposition: %d %d"
            % (len(water.full_system) // water.num_sites, total_water_and_salt)
        )

    # translate so that the water box is centered at the origin
    water.full_system.centerAtOrigin()
    water_residues = water.full_system.splitByResidue()


    for w in water_residues:
        w.reimage()

    sys.stderr.write("Beginning bump-checking water against protein\n")
    # bump-check the water against the protein
    # First step is selecting water and protein heavy atoms
    # These selections can't fail, so we won't trap for zero selection
    # size.
    water_oxygens = loos.selectAtoms(water.full_system, 'name =~ "^O"')
    protein_heavy = loos.selectAtoms(system, '!(name =~ "^H")')

    # find water oxygens within 1.75 Ang of any lipid heavy atom
    clashing_oxygens = water_oxygens.within(1.75, protein_heavy, config.box)
    sys.stderr.write("Found %d clashing waters\n" % (len(clashing_oxygens)))

    # loop over the clashing oxygens, and find which residue each is in,
    # and remove the corresponding residues from the full water box
    for ox in clashing_oxygens:
        i = 0
        for w in water_residues:
            if ox in w:
                water.full_system.remove(w)
                break
        i += 1

    # verify we have enough water
    if len(water.full_system) // water.num_sites < total_water_and_salt:
        raise ValueError(
            "Too few waters after superposition: %d %d"
            % (len(water.full_system) // water.num_sites, total_water_and_salt)
        )

    sys.stderr.write("Finished bump-checking water against protein\n")
    sys.stderr.write(
        "Current # water molecules: %d\n" % (len(water.full_system) // water.num_sites)
    )
    sys.stderr.write("Adding salt\n")

    # regenerate the list of oxygens
    water_oxygens = loos.selectAtoms(water.full_system, 'name =~ "^O"')
    water_residues = water.full_system.splitByResidue()

    # now replace waters with salt
    salts = []
    for salt in config.salt:
        ions = loos.AtomicGroup()
        for i in range(salt.numres):
            a = loos.Atom()
            a.resname(salt.resname)
            a.segid(salt.segname)
            a.name(salt.atomname)
            a.resid(i + 1)

            # pick a water oxygen at random, replace it with salt
            ox = random.choice(water_oxygens)
            a.coords(ox.coords())
            ions.append(a)

            # remove this water
            for w in water_residues:
                if ox in w:
                    water.full_system.remove(w)
                    water_oxygens.remove(ox)
                    break

        salts.append(ions)

    # verify we have enough water
    num_waters = len(water.full_system) // water.num_sites
    if num_waters < config.water.numres:
        raise ValueError(
            "Too few waters after exchanging salt: %d %d"
            % (num_waters, config.water.numres)
        )

    # if we have too many waters, need to delete the difference
    if num_waters > config.water.numres:
        water_residues = water.full_system.splitByResidue()
        diff = num_waters - config.water.numres
        sys.stderr.write("Removing %d waters to get to target size\n" % diff)
        for i in range(diff):
            # remove this water
            ox = random.choice(water_oxygens)
            for w in water_residues:
                if ox in w:
                    water.full_system.remove(w)
                    water_oxygens.remove(ox)
                    break

    # renumber the residues
    for i in range(len(water.full_system)):
        res = i // water.num_sites + 1
        water.full_system[i].resid(res)

    # Replace some of the waters with the internal waters from the protein.
    if config.protein is not None and config.protein.has_water:
        if config.water.numres < len(config.protein.water_seg()) // water.num_sites:
            raise ValueError(
                "Protein has more internal waters than the total target: %d %d"
                % (config.water.numres, len(config.protein.water_seg()) // water.num_sites)
            )
        water.full_system.copyCoordinatesFrom(
            config.protein.water_seg(), 0, len(config.protein.water_seg())
        )


    sys.stderr.write("Assembling final system and writing out coordinates\n")
    # append water and salt to the full system
    system.append(water.full_system)
    for s in salts:
        system.append(s)

    system.renumber()
    system.periodicBox(config.box)

    # Note: all of the manipulation we do mangles the bond list.  Since
    #       we've got a psf anyway, we can safely remove the bond
    #       information here.
    system.clearBonds()

    # write out a final pdb file
    final_pdbfile = open(os.path.join(config.directory, "final.pdb"), "w")
    final_pdb = loos.PDB.fromAtomicGroup(system)
    final_pdbfile.write(str(final_pdb))
    final_pdbfile.close()

    # while we're at it, write out a script to generate a full psf as
    # well, then run psfgen  for them
    psf_script = open(os.path.join(config.directory, "generate_full_psf.inp"), "w")

    psfgen_script = config.generate_psf(True, True, True, False)
    psf_script.write(psfgen_script)
    psf_script.close()


    # If there's a protein, the above psfgen script will depend on the protein's
    # psf file, so we should copy that into the directory as well
    if config.protein is not None:
        shutil.copy(config.protein.psf_file, config.directory)

    # Now, run the psfgen script from the output directory
    sys.stderr.write("Running psfgen to generate final psf\n")
    os.chdir(config.directory)
    psfgen = PSFGen.PSFGen(psfgen_script, config.psfgen_binary)
    psfgen.run()


    # now, read the psf back in and check the total charge
    full_system = loos.createSystem(config.psfname)
    total_charge = full_system.totalCharge()
    if abs(total_charge) > 1e-3:
        sys.stderr.write("\nWARNING WARNING WARNING WARNING\n")
        sys.stderr.write("System has a net charge of %f\n" % total_charge)
        sys.stderr.write(
            "This will likely cause pressure artifacts when you try to run with Ewald.\n"
        )

if __name__ == "__main__":
    main()