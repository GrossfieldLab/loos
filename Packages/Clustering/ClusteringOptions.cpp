#include "ClusteringOptions.hpp"
#include "ClusteringUtils.hpp"
#include <iostream>
#include <string>

namespace Clustering{
  void ClusteringOptions::addGeneric(po::options_description & opts){
    opts.add_options()
    ("similarities,S", po::value<std::string>(&similarity_filename), "File containing whitespace-delimited pairwise similarities.");
  }

  bool ClusteringOptions::postConditions(po::variables_map& vm) {
  	if (similarity_filename.empty())
  		similarityScores = readMatrixFromStream<dtype>(std::cin);
    else {
      std::ifstream stream;
      stream.open(similarity_filename);
      similarityScores = readMatrixFromStream<dtype>(stream);
    }
    return true;
  }

}