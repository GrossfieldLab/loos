add_executable(hcorrelation hcorrelation.cpp hcore.cpp)
target_link_libraries(hcorrelation loos)


add_executable(hbonds hbonds.cpp hcore.cpp)
target_link_libraries(hbonds loos)


add_executable(hmatrix hmatrix.cpp hcore.cpp)
target_link_libraries(hmatrix loos)


add_executable(hcontacts hcontacts.cpp hcore.cpp)
target_link_libraries(hcontacts loos)

