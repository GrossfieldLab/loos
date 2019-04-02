// nmrclust.cpp follows Kelly, Gardner, and Sutcliffe, Prot. Eng. 9 11 1063-1065 (1996) 
// hereafter KGS
// To perform exactly the analysis specified there, one must first apply
// one of the all-to-all rmsd tolls (such as rmsds or multi-rmsds) before
// running this script. Those tools write their RMSD matricies to cout, and
// this script reads from cin, so the effect can be achieved through a pipe.
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "cluster.hpp"
using namespace Eigen;
using namespace std;


int main()
{
  MatrixXd similarity_scores;
  // read similarity scores into matrix, this is step 1 in CGS fig. 1
  similarity_scores = readMatrixFromStream(cin);

  // perform average linkage clustering, step 2 CGS f1.
  // use vector of vectors of structure indices as cluster list.
  vector<vector<uint>> cluster_list;
  for(uint i = 0; i < similarity_scores.rows(); i++){
    vector <uint> cluster{i};
    cluster_list.push_back(cluster);
  }
  
}