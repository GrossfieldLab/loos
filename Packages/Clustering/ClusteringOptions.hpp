#if !defined(LOOS_CLUSTERING_OPTIONS)
#define LOOS_CLUSTERING_OPTIONS

#include "ClusteringTypedefs.hpp"
// #include <loos.hpp>
#include <OptionsFramework.hpp>

namespace Clustering {

namespace opts = loos::OptionsFramework;
namespace po = boost::program_options;
class ClusteringOptions : public opts::OptionsPackage {
public:
  ClusteringOptions() : similarity_filename(""), stream_mode(true) {}
  ClusteringOptions(std::string &similarityFN)
      : similarity_filename(similarityFN), stream_mode(false) {}

  std::string similarity_filename;
  bool stream_mode;
  Eigen::Matrix<dtype, Eigen::Dynamic, Eigen::Dynamic> similarityScores;

private:
  void addGeneric(po::options_description &opts);
  bool postConditions(po::variables_map &map);
};
} // namespace Clustering
#endif