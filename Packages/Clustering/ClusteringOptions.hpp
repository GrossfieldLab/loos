#if !defined(LOOS_CLUSTERING_OPTIONS)
#define LOOS_CLUSTERING_OPTIONS

#include "OptionsFramework.hpp"
#include "ClusteringTypedefs.hpp"

namespace Clustering {

class ClusteringOptions : public OptionsPackage {
public:
  ClusteringOptions(): similarity_filename("") {}
  ClusteringOptions(std::string& similarityFN): similarity_filename(similarityFN) {}

  std::string similarity_filename;
  Eigen::Matrix<dtype, Eigen::Dynamic, Eigen::Dynamic> similarityScores;



private:
  void addGeneric(po::options_description &opts);
  bool postConditions(po::variables_map& map);

};

} // namespace Clustering
#endif