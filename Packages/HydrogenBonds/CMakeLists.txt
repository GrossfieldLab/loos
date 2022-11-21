add_executable(hcorrelation hcorrelation.cpp hcore.cpp)
target_link_libraries(hcorrelation loos)
install(TARGETS hcorrelation)


add_executable(hbonds hbonds.cpp hcore.cpp)
target_link_libraries(hbonds loos)
install(TARGETS hbonds)


add_executable(hmatrix hmatrix.cpp hcore.cpp)
target_link_libraries(hmatrix loos)
install(TARGETS hmatrix)


add_executable(hcontacts hcontacts.cpp hcore.cpp)
target_link_libraries(hcontacts loos)
install(TARGETS hcontacts)

