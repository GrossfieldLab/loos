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
  cout << similarity_scores << '\n';
  removeRow(similarity_scores, 2);
  removeCol(similarity_scores, 1);
  cout << similarity_scores << '\n';
  
}