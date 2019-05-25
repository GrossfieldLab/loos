// cluster-kgs.cpp 
// Kelly, Gardner, and Sutcliffe, Prot. Eng. 9 11 1063-1065 (1996) 
// hereafter KGS
// To perform exactly the analysis specified there, one must first apply
// one of the all-to-all rmsd tolls (such as rmsds or multi-rmsds) before
// running this tool. Those tools write their RMSD matricies to cout, and
// this tool reads from cin, so the effect can be achieved through a pipe.
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "Clustering.hpp"

using namespace Eigen;
using namespace std;
using namespace Clustering;

std::string helpstr = "XXX";

int main(int argc, char* argv)
{
  if (argc < 2)
  {
    cout << helpstr << endl;
    exit(0);
  }
  else
  {
    for (uint i = 0; i < argc; i++)
    {
      if (argv[i] == "-h" || argv[i] == "--help")
      {
        cout << helpstr << endl;
        exit(0);
      }
    }
  }
  MatrixXd similarityScores = readMatrixFromStream(cin);
  KGS clusterer(similarityScores);
  clusterer.cluster();
  uint optStg = clusterer.cutoff();
  vector<uint> exemplars = getExemplars(clusterer.clusterTraj[optStg], clusterer.refDists);
  // below here is output stuff. All quantities of interest have been obtained.
  string indent = "  ";
  string offset = "    ";
  cout << "{";
  cout << indent + "\"optimal stage\": "<< optStg << "," << endl;
  cout << indent + "\"penalties\": ";
  containerAsInlineJSONArr(clusterer.penalties, cout);
  cout << "," << endl;
  cout << indent + "\"clusters\": ";
  vectorVectorAsJSONArr(clusterer.clusterTraj[optStg], cout, offset=offset);
  cout << "," << endl;
  cout << indent + "\"exemplars\": ";
  containerAsJSONArr(exemplars, cout, offset=offset);
  cout << endl;
  cout << "}";
}