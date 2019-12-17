#include "ClusteringOptions.hpp"

namespace Clustering{
  void ClusteringOptions::addGeneric(po::options_description & opts){
    opts.add_options()
    ("similarities,S", po::value<std::string>(&similarity_filename), "File containing whitespace-delimited pairwise similarities.");
  }

  // postConditions might not be needed yet for such a crude file.

}