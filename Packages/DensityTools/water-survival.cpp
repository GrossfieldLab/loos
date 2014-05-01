// (c) 2014 Tod D. Romo, Grossfield Lab, URMC

#include <loos.hpp>

using namespace loos;
using namespace std;

typedef Math::Matrix<int>   WaterClassMatrix;



int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  if (argc == 1) {
    cerr << "Usage- " << argv[0] << " water_matrix [max-t] >output.asc\n";
    exit(-1);
  }

  int k = 1;
  string matname(argv[k++]);
  uint max_t = 0;
  if (k != argc)
    max_t = strtoul(argv[k++], 0, 10);
  
  Math::Matrix<int> M;
  cerr << "Reading matrix...\n";
  readAsciiMatrix(matname, M);
  uint m = M.rows();
  uint n = M.cols();

  if (max_t == 0)
    max_t = n/10;

  cerr << boost::format("Water matrix is %d x %d\n") % m % n;

  cout << "# " << hdr << endl;
  cout << "# tau\tavg\tstdev\tsterr\n";
  
  cerr << "Processing- ";
  

  for (uint tau=1; tau<max_t; ++tau) {
    vector<double> survivals;
    if (tau % 100 == 0)
      cerr << '.';
    
    for (uint j=0; j<m; ++j) {
      uint inside = 0;
      uint pairs = 0;
      
      for (uint t=0; t<n-tau-1; ++t)
	if (M(j, t))
	{
	  ++pairs;
	  if (M(j, t+tau))
	    ++inside;
	}
      if (pairs)
	survivals.push_back(static_cast<double>(inside) / (pairs));
    }
    dTimeSeries ts(survivals);
    cout << tau << '\t' << ts.average() << '\t' << ts.stdev() << '\t' << ts.sterr() << endl;
  }
  
  cerr << " Done\n";

}



  
