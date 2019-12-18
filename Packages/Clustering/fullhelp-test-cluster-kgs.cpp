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
  "usage: \n"
"cluster-kgs < similarity_scores.asc > clustering_results.json \n"
" \n"
"kgsclus mimics the clustering aspect of the NMRCLUST utility that is \n"
"incorporated as part of the OLDERADO webserver for structural biology \n"
"informatics. It was originally published as: \n"
"Kelly, Gardner, and Sutcliffe, Prot. Eng. 9 11 1063-1065 (1996) \n"
"This type of clustering exists in other places, most notably in R, and has been\n"
" put to many uses beside clustering protein structures with their RMSD as the \n"
"distance between each structure. It is called cluster-kgs because this method \n"
"is referred to in other contexts (that is, where it is not being used to \n"
"analyze NMR ensembles) as kgs clustering, and because this executable operates \n"
"on a provided similarity matrix it is similarly flexible. Note that we do not \n"
"implement the 'eigen analysis' for cluster center determination, instead \n"
"choosing to use the element from each cluster with the lowest mean distance to \n"
"the other elements in the cluster.  \n"
" \n"
"The tool works by reading in a similarity score matrix from a file (or stdin) \n"
"and writing the clustering results to stdout. The results report the index of \n"
"each cluster, with all the elements in each cluster following its index on the \n"
"same line. It will also provide an exemplar (the element nearest the centroid) \n"
"for each cluster in a separate block. The input matrix should be an NxN \n"
"symmetric matrix of similarity scores where the ij-th element is the similarity\n"
" between datum i and datum j. The similarity score matrix is expected to be \n"
"whitespace delimited in the column and newline delimited in the row. '#' is an \n"
"acceptable comment character, but only produces a comment-read at the beginning\n"
" of a line (not at any point in a line, like a comment in a shell script). \n"
" \n"
"In order to mimic the functionality of the OLDERADO tool mentioned above, one \n"
"can use the loos tool rmsds (or similar) to produce the matrix of similarity \n"
"scores.  \n"
"For example:  \n"
" \n"
"$ rmsds model.pdb ensemble.dcd | cluster-kgs -s > clustering_results.json \n"
" \n"
"would use rmsds to compute the alpha carbon RMSDs from the frame-pairs in \n"
"ensemble.dcd to generate the similarity matrix, then redirect it to cluster-\n"
"kgs, which will read from stdin because the -s flag was thrown. Then the \n"
"clustering results are written to an output file (which should be valid JSON, \n"
"for convenient further scripting). This shell-redirect would also cause the \n"
"distance matrix from rmsds to be written to stderr. Note that in this \n"
"particular command line the RMSD values emitted by rmsds will be in angstroms, \n"
"and will be rounded to 2 digits. For more reported precision (rmsds uses \n"
"doubles internally), use the '-p' flag. If you would like to both save the \n"
"similarities generated in this way, but also not have them written to disk \n"
"before feeding them to the clustering algorithm, you can redirect stderr and \n"
"stdout to separate files: \n"
" \n"
"$ rmsds model.pdb ensemble.dcd | cluster-kgs -s 1> clustering_results.json \n"
"2>distances.asc \n"
" \n"
"You can also read a distance matrix from a file using the -f flag. If you do \n"
"that, it will not be emitted to stderr, and you would write: \n"
"cluster-kgs -f distances.asc > clustering_results.json \n"
" \n"
"The output from the clustering will be structured as JSON, and will have four \n"
"keys. One will be the invocation header, to help record where the data came \n"
"from. Then one recording which stage of the clustering was chosen as the cutoff\n"
" stage. One recording the penalty values at each stage (they could potentially \n"
"be quite similar). One will be a list (in no particular order) of the clusters,\n"
" where each cluster will be a list of element indexes corresponding to the \n"
"members of that cluster. Finally the last will be a list of exemplars for each \n"
"cluster. These will be in the same order as the list of clusters, so the first \n"
"exemplar in the list is the exemplar for the first cluster in the cluster list.\n"
;
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