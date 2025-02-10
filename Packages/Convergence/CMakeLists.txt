set(FIDLIB fid-lib.cpp bcomlib.hpp fid-lib.hpp)

add_executable(bcom bcom.cpp ${FIDLIB})
target_link_libraries(bcom loos)
install(TARGETS bcom)

add_executable(boot_bcom boot_bcom.cpp ${FIDLIB})
target_link_libraries(boot_bcom loos)
install(TARGETS boot_bcom)

add_executable(ufidpick ufidpick.cpp ${FIDLIB})
target_link_libraries(ufidpick loos)
install(TARGETS ufidpick)

add_executable(assign_frames assign_frames.cpp ${FIDLIB})
target_link_libraries(assign_frames loos)
install(TARGETS assign_frames)

add_executable(decorr_time decorr_time.cpp ${FIDLIB})
target_link_libraries(decorr_time loos)
install(TARGETS decorr_time)

add_executable(coscon coscon.cpp ${FIDLIB})
target_link_libraries(coscon loos)
install(TARGETS coscon)

add_executable(qcoscon qcoscon.cpp ${FIDLIB})
target_link_libraries(qcoscon loos)
install(TARGETS qcoscon)

add_executable(rsv-coscon rsv-coscon.cpp ${FIDLIB})
target_link_libraries(rsv-coscon loos)
install(TARGETS rsv-coscon)

add_executable(fidpick fidpick.cpp)
target_link_libraries(fidpick loos)
install(TARGETS fidpick)

add_executable(sortfids sortfids.cpp)
target_link_libraries(sortfids loos)
install(TARGETS sortfids)

add_executable(hierarchy hierarchy.cpp)
target_link_libraries(hierarchy loos)
install(TARGETS hierarchy)

add_executable(neff neff.cpp)
target_link_libraries(neff loos)
install(TARGETS neff)

add_executable(block_average block_average.cpp)
target_link_libraries(block_average loos)
install(TARGETS block_average)

add_executable(avgconv avgconv.cpp)
target_link_libraries(avgconv loos)
install(TARGETS avgconv)

add_executable(block_avgconv block_avgconv.cpp)
target_link_libraries(block_avgconv loos)
install(TARGETS block_avgconv)

add_executable(expfit expfit.cpp)
target_link_libraries(expfit loos)
install(TARGETS expfit)

add_executable(chist chist.cpp)
target_link_libraries(chist loos)
install(TARGETS chist)

### Handle scripts
install(PROGRAMS effsize.pl bootstrap_overlap.pl DESTINATION bin)