// nmrclust.cpp 
// Kelly, Gardner, and Sutcliffe, Prot. Eng. 9 11 1063-1065 (1996) 
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
  cout << "optimal stage:  " << optStg << endl;
  cout << "penalties:  " << clusterer.penalties << endl;
  clusterer.writeClusters(optStg, cout);
  vector<uint> exemplars = getExemplars(clusterer.clusterTraj[optStg], clusterer.refDists);
  // print exemplars out below here
  cout << "Exemplars:  " << endl;
  for (uint i = 0; i < exemplars.size(); i++)
    cout << i << ' ' << exemplars[i] << endl;
}