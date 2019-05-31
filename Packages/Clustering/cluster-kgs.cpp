// cluster-kgs.cpp
// Kelly, Gardner, and Sutcliffe, Prot. Eng. 9 11 1063-1065 (1996)
// hereafter KGS
// To perform exactly the analysis specified there, one must first apply
// one of the all-to-all rmsd tolls (such as rmsds or multi-rmsds) before
// running this tool. Those tools write their RMSD matricies to cout, and
// this tool reads from cin, so the effect can be achieved through a pipe.
#include "Clustering.hpp"
#include <iostream>
#include <string>

using std::cout;
using std::cin;
using std::endl;
using std::vector;

// using namespace Clustering;

const std::string helpstr = "XXX";
const std::string indent = "  ";
int main(int argc, char* argv[])
{
  if (argc < 2) {
    cout << helpstr << endl;
    exit(0);
  } else {
    for (int i = 0; i < argc; i++) {
      if (argv[i] == "-h" || argv[i] == "--help") {
        cout << helpstr << endl;
        exit(0);
      }
    }
  }
  Eigen::MatrixXd similarityScores = Clustering::readMatrixFromStream(cin);
  cout << similarityScores;
  // Clustering::KGS clusterer(similarityScores);
  // clusterer.cluster();
  // uint optStg = clusterer.cutoff();
  // vector<uint> exemplars =
  //   getExemplars(clusterer.clusterTraj[optStg], clusterer.refDists);
  // // below here is output stuff. All quantities of interest have been obtained.
  // cout << "{";
  // cout << indent + "\"optimal stage\": " << optStg << "," << endl;
  // cout << indent + "\"penalties\": ";
  // containerAsOneLineJSONArr<Eigen::VectorXd>(clusterer.penalties, cout);
  // cout << "," << endl;
  // cout << indent + "\"clusters\": ";
  // vectorVectorsAsJSONArr<uint>((clusterer.clusterTraj)[optStg], cout, "  ");
  // cout << "," << endl;
  // cout << indent + "\"exemplars\": ";
  // containerAsJSONArr<vector<uint>>(exemplars, cout, "  ");
  // cout << endl;
  // cout << "}";
}