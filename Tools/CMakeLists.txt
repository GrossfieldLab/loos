find_package(Boost REQUIRED COMPONENTS thread filesystem)

set(LOOS_BASIC_TOOLS
    aligner
    area_per_lipid
    atomic-rdf
    averager
    big-svd
    bounding
    center-molecule
    clipper
    concat-selection
    contact-time
    contacts
    convert2pdb
    coverlap
    cross-dist
    crossing-waters
    dcdinfo
    density-dist-windowed
    density-dist
    dibmops
    dihedrals
    distance_prob
    drifter
    dumpmol
    esp_mesh
    exposure
    fcontacts
    fixdcd
    frame2pdb
    helix_kink
    hex-order
    interdist
    kurskew
    lipid_survival
    long-bond-finder
    loos_version
    membrane_map
    merge-traj
    model-meta-stats
    model-select
    model2matlab
    molshape
    mops
    native_contacts
    ocf
    order_params
    packing_score
    paxes
    periodic_box
    perturb-structure
    phase-pdb
    porcupine
    rad-gyr
    ramachandran
    rdf
    rebond
    recenter-trj
    reimage-by-molecule
    renum-pdb
    residue-contact-map
    rgyr
    rmsd2ref
    rmsf
    rmsfit
    rna_suites
    rotamer
    serialize-selection
    smooth-traj
    subsetter
    svdcolmap
    torsion
    traj2dcd
    traj2matlab
    traj2pdb
    trajinfo
    transition_contacts
    verap
    xtcinfo
    xy_rdf
    xy_rdf_timeseries

)

foreach(TOOL ${LOOS_BASIC_TOOLS})
  add_executable(${TOOL} ${TOOL}.cpp)
  target_link_libraries(${TOOL} loos)
  install(TARGETS ${TOOL})
endforeach()

add_executable(multi-rmsds multi-rmsds.cpp)
target_link_libraries(multi-rmsds loos Boost::thread)
install(TARGETS multi-rmsds)


add_executable(rms-overlap rms-overlap.cpp)
target_link_libraries(rms-overlap loos Boost::thread)
install(TARGETS rms-overlap)

add_executable(rmsds rmsds.cpp)
target_link_libraries(rmsds loos Boost::thread)
install(TARGETS rmsds)

add_executable(svd svd.cpp)
target_link_libraries(svd loos Boost::filesystem)
install(TARGETS svd)


### Handle scripts
install(PROGRAMS blurrogram.pl  gmxdump2pdb.pl mat2loos.pl DESTINATION bin)