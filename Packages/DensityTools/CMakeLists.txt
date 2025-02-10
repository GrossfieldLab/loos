add_library(loos_density
	    GridUtils.cpp internal-water-filter.cpp water-hist-lib.cpp water-lib.cpp
	    DensityGrid.hpp GridUtils.hpp internal-water-filter.hpp water-hist-lib.hpp water-lib.hpp DensityOptions.hpp
	    )

target_include_directories(loos_density PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(loos_density loos)
install(TARGETS loos_density DESTINATION lib)

add_executable(gridinfo gridinfo.cpp)
target_link_libraries(gridinfo loos_density)
install(TARGETS gridinfo)

add_executable(grid2ascii grid2ascii.cpp)
target_link_libraries(grid2ascii loos_density)
install(TARGETS grid2ascii)

add_executable(grid2xplor grid2xplor.cpp)
target_link_libraries(grid2xplor loos_density)
install(TARGETS grid2xplor)

add_executable(gridgauss gridgauss.cpp)
target_link_libraries(gridgauss loos_density)
install(TARGETS gridgauss)

add_executable(gridscale gridscale.cpp)
target_link_libraries(gridscale loos_density)
install(TARGETS gridscale)

add_executable(gridslice gridslice.cpp)
target_link_libraries(gridslice loos_density)
install(TARGETS gridslice)

add_executable(gridmask gridmask.cpp)
target_link_libraries(gridmask loos_density)
install(TARGETS gridmask)

add_executable(blobid blobid.cpp)
target_link_libraries(blobid loos_density)
install(TARGETS blobid)

add_executable(contained contained.cpp)
target_link_libraries(contained loos_density)
install(TARGETS contained)

add_executable(gridstat gridstat.cpp)
target_link_libraries(gridstat loos_density)
install(TARGETS gridstat)

add_executable(peakify peakify.cpp)
target_link_libraries(peakify loos_density)
install(TARGETS peakify)

add_executable(pick_blob pick_blob.cpp)
target_link_libraries(pick_blob loos_density)
install(TARGETS pick_blob)

add_executable(blob_stats blob_stats.cpp)
target_link_libraries(blob_stats loos_density)
install(TARGETS blob_stats)

add_executable(water-inside water-inside.cpp)
target_link_libraries(water-inside loos_density)
install(TARGETS water-inside)

add_executable(water-extract water-extract.cpp)
target_link_libraries(water-extract loos_density)
install(TARGETS water-extract)

add_executable(water-hist water-hist.cpp)
target_link_libraries(water-hist loos_density)
install(TARGETS water-hist)

add_executable(water-count water-count.cpp)
target_link_libraries(water-count loos_density)
install(TARGETS water-count)

add_executable(water-sides water-sides.cpp)
target_link_libraries(water-sides loos_density)
install(TARGETS water-sides)

add_executable(blob_contact blob_contact.cpp)
target_link_libraries(blob_contact loos_density)
install(TARGETS blob_contact)

add_executable(griddiff griddiff.cpp)
target_link_libraries(griddiff loos_density)
install(TARGETS griddiff)

add_executable(near_blobs near_blobs.cpp)
target_link_libraries(near_blobs loos_density)
install(TARGETS near_blobs)

add_executable(gridautoscale gridautoscale.cpp)
target_link_libraries(gridautoscale loos_density)
install(TARGETS gridautoscale)

add_executable(gridavg gridavg.cpp)
target_link_libraries(gridavg loos_density)
install(TARGETS gridavg)

add_executable(water-autocorrel water-autocorrel.cpp)
target_link_libraries(water-autocorrel loos_density)
install(TARGETS water-autocorrel)

add_executable(water-survival water-survival.cpp)
target_link_libraries(water-survival loos_density)
install(TARGETS water-survival)

### Handle scripts
install(PROGRAMS grid_slicer.pl DESTINATION bin)