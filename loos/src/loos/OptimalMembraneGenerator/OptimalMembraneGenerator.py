#!/usr/bin/env python3

import random

import loos
import loos.OptimalMembraneGenerator
from loos.OptimalMembraneGenerator import PSFGen
from loos.OptimalMembraneGenerator import NAMD
from loos.OptimalMembraneGenerator import WaterBox

def main():
    import sys
    import shutil
    import os
    import stat

    config = PSFGen.ReadConfig(sys.argv[1])
    if len(sys.argv) > 2:
        config.directory = sys.argv[2]

    box_scaling = 3.0
    box = box_scaling * config.box
    # ridiculously large z dimension used initially, chosen to be effectively infinite
    huge_z = 10000.0

    box.z(huge_z)
    scaled_box = box.x()
    half_box = 0.5 * scaled_box

    # the distance that is considered a "contact"
    overlap_dist = 3.0
    # the number of contacts required for a molecule to be rejected
    overlap_threshold = 3

    # set how big the "water" region of the box will be
    water_width = config.box.z()

    # Don't need a blank spot in the middle if there are no lipids
    if len(config.segments) > 0:
        water_width -= 29.0  # TODO: need a smarter algorithm
        #       for this, probably using
        #       the phos values

    # set the number of minimization and dynamics NAMD will do on each cycle
    number_of_steps = 100

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

    # If the "protein" is actually a bunch of independent molecules
    # (e.g. a bunch of peptides), we'll want to scale them in x & y
    # to match the expanded box.
    if config.protein is not None and config.protein.scale:
        protein_molecules = config.protein.model.splitByMolecule()
        for m in protein_molecules:
            centroid = m.centroid()
            scaled = box_scaling * centroid
            diff = scaled - centroid
            diff.z(0.0)  # we only want to translate in the xy plane
            m.translate(diff)

    molecules = system.splitByMolecule()

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

    x_axis = loos.GCoord(1, 0, 0)
    z_axis = loos.GCoord(0, 0, 1)

    for j in range(len(segments)):
        seg_ag = segments[j]
        seg_ag_arr = seg_ag.splitByMolecule()
        seg_conf = config.segments[j]
        i = 0

        while i < seg_conf.numres:
            lipid = seg_conf.library.pick_structure()
            phos = loos.selectAtoms(lipid, 'name == "' + seg_conf.phos_atom + '"')
            if len(phos) == 0:
                sys.stderr.write(
                    'Selection failed: "phos" atom %s doesn\'t exist in lipid %s\n'
                    % (seg_conf.phos_atom, seg_conf.resname)
                )
                sys.stderr.write("Exiting...\n")
                sys.exit(0)

            # put the molecule at the origin
            lipid.centerAtOrigin()

            # perform random rotation about z axis
            lipid.rotate(z_axis, random.uniform(0.0, 360.0))

            if seg_conf.placement < 0:
                lipid.rotate(x_axis, 180)

            # put the phosphate back at the origin
            lipid.translate(-phos[0].coords())

            # generate the new coordinates
            x = (random.random() * scaled_box) - half_box
            y = (random.random() * scaled_box) - half_box

            # set the z location, with +/- 0.5 random noise
            z = seg_conf.phos_height + (random.random() - 0.5)

            if seg_conf.placement < 0:
                z = -z

            vec = loos.GCoord(x, y, z)
            lipid.translate(vec)

            # do a bump-check against previously placed lipids
            # If at any time we realize we've got too many clashes,
            # we'll break out of the loop early, and try to place
            # this lipid again by bouncing to the top of the while
            # loop without incrementing the index i

            # first, bump check against the protein, if any
            num_overlap = 0
            if config.protein is not None:
                for seg in config.protein.segments:
                    overlap = lipid.within(overlap_dist, seg, box)
                    num_overlap += overlap.size()
                    if num_overlap > overlap_threshold:
                        break
                if num_overlap > overlap_threshold:
                    continue

            # second, check this lipid against all previous segments
            for k in range(j):
                overlap = lipid.within(overlap_dist, segments[k], box)
                num_overlap += overlap.size()
                if num_overlap > overlap_threshold:
                    break
            if num_overlap > overlap_threshold:
                continue

            # now, check earlier lipids in this segment
            for k in range(i):
                overlap = lipid.within(overlap_dist, seg_ag_arr[k], box)
                num_overlap += overlap.size()
                if num_overlap > overlap_threshold:
                    break
            if num_overlap > overlap_threshold:
                continue

            # copy the coordinates back into the real object
            seg_ag_arr[i].copyMappedCoordinatesFrom(lipid)

            i += 1

    # copy the protein coordinates into the system
    if config.protein is not None:
        for s in config.protein.segments:
            current_seg = s[0].segid()
            # Don't need to trap failed selection here, because we
            # already know this segment exists
            seg = loos.selectAtoms(system, 'segname == "' + current_seg + '"')
            seg.copyMappedCoordinatesFrom(s)

    # TODO: decide if we want to put this back. Assume
    #       the "protein" is already centered"
    # system.centerAtOrigin()
    pdb = loos.PDB.fromAtomicGroup(system)
    pdb_out = open(os.path.join(config.directory, "lipid_only.pdb"), "w")
    pdb_out.write(str(pdb))
    pdb_out.close()

    # Now start scaling the box down with successive minimizations
    num_iter = 100
    delta = (scaled_box - config.box.x()) / num_iter

    # only need to do the scaling and minimization if we've got lipids
    if len(config.segments) > 0:

        next_pdb = "lipid_only.pdb"
        dim = scaled_box
        for i in range(num_iter + 1):
            old_dim = dim
            dim = scaled_box - (i * delta)

            # scale the centers of mass of individual molecules
            # TODO: this will have to change to handle "fixed" molecules
            scale_factor = dim / old_dim
            for j in range(len(molecules)):
                centroid = molecules[j].centroid()
                centroid.z(0.0)  # only scale in plane of membrane
                new_centroid = centroid * scale_factor
                diff = new_centroid - centroid
                molecules[j].translate(diff)

            # output the new coordinates so we can minimize with NAMD
            pdb_out = open(os.path.join(config.directory, next_pdb), "w")
            pdb_out.write(str(pdb))
            pdb_out.close()

            # TODO: long term, I'd like to change this so that instead of
            #       running a new instance of NAMD, we're reusing the same
            #       process, and adding new commands
            current_box = loos.GCoord(dim, dim, huge_z)
            core_name = "lipid_shrink_" + str(i)
            full_core = os.path.join(config.directory, core_name)

            # copy the psf file to the working directory
            local_psfname = os.path.basename(config.psfname)

            namd = NAMD.NAMD(
                local_psfname,
                next_pdb,
                core_name,
                config.parameters,
                current_box,
                config.namd_binary,
            )
            namd_inputfilename = full_core + ".inp"
            namd_outputfilename = full_core + ".out"
            namd.write_inputfile(namd_inputfilename, number_of_steps)

            namd.write_restraintfile(config.directory, system)
            namd.run_namd(namd_inputfilename, namd_outputfilename)

            # copy the final coordinates to a pdb file
            shutil.copyfile(
                full_core + ".coor", os.path.join(config.directory, next_pdb)
            )

            # read in the new structure
            new_system = loos.createSystem(os.path.join(config.directory, next_pdb))

            # copy coordinates back into old AtomicGroup
            system.copyCoordinatesFrom(new_system, 0, system.size())

            # TODO: this will have to change to handle "fixed" molecules
            system.centerAtOrigin()
            sys.stderr.write("lipid only minimization cycle: %d\n" % i)

    pdb_out = open(os.path.join(config.directory, "lipid_min.pdb"), "w")
    pdb_out.write(str(pdb))
    pdb_out.close()

    sys.stderr.write("Beginning water box construction\n")
    # now add water and salt
    water_template = loos.GCoord(1.0, 1.0, 1.0)
    water_template *= config.water.box_size
    water_target = loos.GCoord(config.box.x(), config.box.y(), water_width)
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

    # translate so that the water box is centered on the periodic boundary,
    # then set its periodic box to the full system box size and reimage
    water.full_system.periodicBox(config.box)
    trans = loos.GCoord(0.0, 0.0, 0.5 * config.box.z())
    water.full_system.translate(trans)
    water_residues = water.full_system.splitByResidue()

    for w in water_residues:
        w.reimage()

    sys.stderr.write("Beginning bump-checking water against lipid\n")
    # bump-check the water against the lipids
    # First step is selecting water and lipid heavy atoms
    # These selections can't fail, so we won't trap for zero selection
    # size.
    water_oxygens = loos.selectAtoms(water.full_system, 'name =~ "^O"')
    lipid_heavy = loos.selectAtoms(system, '!(name =~ "^H")')

    # find water oxygens within 1.75 Ang of any lipid heavy atom
    clashing_oxygens = water_oxygens.within(1.75, lipid_heavy, config.box)
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

    sys.stderr.write("Finished bump-checking water against lipid\n")
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
                % (
                    config.water.numres,
                    len(config.protein.water_seg()) // water.num_sites,
                )
            )
        water.full_system.copyCoordinatesFrom(
            config.protein.water_seg(), 0, len(config.protein.water_seg())
        )
        internal_water = loos.selectAtoms(
            system, 'segid == "' + config.protein.water_segname + '"'
        )
        system.remove(internal_water)

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
            "This will likely cause pressure and area artifacts when you try to run with Ewald.\n"
        )


if __name__ == "__main__":
    main()