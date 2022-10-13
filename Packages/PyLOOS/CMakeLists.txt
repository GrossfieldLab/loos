set(PYLOOS_PACKAGE_SCRIPTS
    add_bonds.py
    all_contacts.py
    all_stacking.py
    axis_with_membrane.py
    center_molecule.py
    cluster-structures.py
    contact_distance.py
    cylindrical-density.py
    cylindrical-thickness.py
    domain.py
    flip_chirality.py
    helix_axes.py
    hierarchical-cluster.py
    inside_helices.py
    native-hbs.py
    packing_score_per_res.py
    protein_tilt.py
    rare-event-detection.py
    renum-by-mol.py
    rmsd_to_average.py
    scattering.py
    set_segids.py
    simple_model_calc.py
    simple_model_transform.py
    simple_traj_calc.py
    simple_traj_transform.py
    sysinfo.py
    total_charge.py
)

#install(PROGRAMS ${PYLOOS_PACKAGE_SCRIPTS} DESTINATION bin)

foreach(script_name  ${PYLOOS_PACKAGE_SCRIPTS})
    string(REGEX REPLACE "^(.+)\.py$" "\\1" program_name ${script_name})
    install(PROGRAMS ${script_name} DESTINATION bin RENAME ${program_name})
endforeach()
