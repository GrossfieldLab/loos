#include "ClusteringOptions.hpp"
#include "ClusteringUtils.hpp"
#include <iostream>
#include <string>

namespace Clustering{
  void ClusteringOptions::addGeneric(po::options_description & opts){
    opts.add_options()
    ("score-file,f", po::value<std::string>(&similarity_filename), "File containing whitespace-delimited pairwise similarities.")
    ("stream,s", po::bool_switch(&stream_mode), "Read similarities from stdin.");
  }

  bool ClusteringOptions::postConditions(po::variables_map& vm) {
    if (stream_mode && !similarity_filename.empty()){
      std::cerr << "Usage Error: \nYou've asked that the tool both reads a file: \n  \"" << similarity_filename << "\"\nand that it read from stdin. \n";
      return false;
    }
      
  	if (stream_mode)
  		similarityScores = readMatrixFromStream<dtype>(std::cin);
    else if (!similarity_filename.empty()){
      std::ifstream stream;
      stream.open(similarity_filename);
      similarityScores = readMatrixFromStream<dtype>(stream);
    }
    else // false causes briefhelp print and exit
      return false;
    
    return true;
  }

}