add_library(loos_enm
	    spring_functions.cpp enm-lib.cpp vsa-lib.cpp
	    anm-lib.hpp enm-lib.hpp spring_functions.hpp vsa-lib.hpp
	    )

target_include_directories(loos_enm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(loos_enm loos)
install(TARGETS loos_enm DESTINATION lib)

add_executable(gnm gnm.cpp)
target_link_libraries(gnm loos)
install(TARGETS gnm)

add_executable(enmovie enmovie.cpp)
target_link_libraries(enmovie loos)
install(TARGETS enmovie)

add_executable(psf-masses psf-masses.cpp)
target_link_libraries(psf-masses loos)
install(TARGETS psf-masses)

add_executable(heavy-ca heavy-ca.cpp)
target_link_libraries(heavy-ca loos)
install(TARGETS heavy-ca)

add_executable(eigenflucc eigenflucc.cpp)
target_link_libraries(eigenflucc loos)
install(TARGETS eigenflucc)


add_executable(vsa vsa.cpp)
target_link_libraries(vsa loos_enm)
install(TARGETS vsa)

add_executable(anm anm.cpp)
target_link_libraries(anm loos_enm)
install(TARGETS anm)
