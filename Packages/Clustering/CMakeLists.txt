find_package(Eigen3 REQUIRED NO_MODULE)

add_executable(cluster-kgs cluster-kgs.cpp
  ClusteringUtils.cpp HAC.cpp AverageLinkage.cpp KGS.cpp ClusteringOptions.cpp
  Clustering.hpp ClusteringUtils.hpp HAC.hpp AverageLinkage.hpp KGS.hpp ClusteringOptions.hpp)

target_link_libraries(cluster-kgs loos Eigen3::Eigen)
install(TARGETS cluster-kgs)

### Handle scripts
install(PROGRAMS cluster_pops.py frame-picker.py DESTINATION bin)