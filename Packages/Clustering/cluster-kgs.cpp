// cluster-kgs.cpp
// Kelly, Gardner, and Sutcliffe, Prot. Eng. 9 11 1063-1065 (1996)
// hereafter KGS
// To perform exactly the analysis specified there, one must first apply
// one of the all-to-all rmsd tolls (such as rmsds or multi-rmsds) before
// running this tool. Those tools write their RMSD matricies to cout, and
// this tool reads from cin, so the effect can be achieved through a pipe.
#include "Clustering.hpp"
#include "ClusteringOptions.hpp"
#include <iostream>
#include <loos.hpp>
#include <string>

using std::cin;
using std::cout;
using std::endl;
using std::vector;

using namespace Clustering;
// namespace opts = loos::OptionsFramework;

const std::string indent = "  ";

std::string fullHelpMessage(void) {
  std::string helpstr = 
  "XXX";
  return helpstr;
}
int main(int argc, char *argv[]) {
  std::string hdr = loos::invocationHeader(argc, argv);
  opts::BasicOptions *bopts = new opts::BasicOptions(fullHelpMessage());
  ClusteringOptions *copts = new ClusteringOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(copts);
  if (!options.parse(argc, argv))
    exit(-1);

  KGS clusterer(copts->similarityScores);
  if (copts->stream_mode)
    std::cerr << copts->similarityScores; // Eigen objects stringify.
  clusterer.cluster();
  idxT optStg = clusterer.cutoff();
  vector<idxT> exemplars =
      getExemplars(clusterer.clusterTraj[optStg], copts->similarityScores);

  // below here is output stuff. All quantities of interest have been obtained.
  cout << "{\n";
  cout << indent + "\"invocation\": \"" << hdr << "\",\n";
  cout << indent + "\"optimal stage\": " << optStg << ",\n";
  cout << indent + "\"penalties\": ";
  containerAsJSONArr<Eigen::Matrix<dtype, Eigen::Dynamic, 1>>(
      clusterer.penalties, cout);
  cout << ",\n";
  cout << indent + "\"clusters\": ";
  vectorVectorsAsJSONArr<idxT>((clusterer.clusterTraj)[optStg], cout, "  ");
  cout << ",\n";
  cout << indent + "\"exemplars\": ";
  containerAsJSONArr<vector<idxT>>(exemplars, cout, "  ");
  cout << "\n";
  cout << "}" << endl;
}