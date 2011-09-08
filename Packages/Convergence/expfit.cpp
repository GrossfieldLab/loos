#include <loos.hpp>
#include <Simplex.hpp>

using namespace std;
using namespace loos;


typedef pair<double, double>     Point;
typedef vector<double>           vecDouble;
typedef vector<vecDouble>        vvecDouble;

struct ExponentialFit {
  ExponentialFit(const vector<Point>& datapoints) : _datapoints(datapoints) { }

  double operator()(const vector<double>& v) {
    double sum = 0.0;

    for (vector<Point>::iterator j = _datapoints.begin(); j != _datapoints.end(); ++j) {
      double x = j->first;
      double y = j->second;

      double val = 1.0;
      for (vector<double>::const_iterator i = v.begin(); i != v.end();) {
        double k = *i++;
        double t = *i++;
        val += k * exp(-x/t);
      }
      double d = y - val;
      sum += d*d;
    }

    return(sum);
  }


  vector<Point> _datapoints;
};


vvecDouble readData(const string& fname) {
  ifstream ifs(fname.c_str());
  return(readTable<double>(ifs));
}



int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);



  int k = 1;
  vvecDouble bcom = readData(argv[k++]);
  vvecDouble bbcom = readData(argv[k++]);

  if (bcom.size() != bbcom.size()) {
    cerr << boost::format("Error- bcom has %d datapoints but bbcom has %d\n") % bcom.size() % bbcom.size();
    exit(-10);
  }
  
  vector<Point> datapoints;
  for (uint i=0; i<bcom.size(); ++i)
    datapoints.push_back(Point(bcom[i][0], bbcom[i][1]/bcom[i][1]));

  int nreps = atoi(argv[k++]);
  int ndims = atoi(argv[k++]);
  vector<double> seeds;
  vector<double> lens;
  ndims *= 2;
  for (int i=0; i<ndims; ++i) {
    double d = strtod(argv[k++], 0);
    seeds.push_back(d);
    lens.push_back(d/2.0);
  }

  ExponentialFit fitFunction(datapoints);

  while (nreps-- > 0) {
    Simplex<double> optimizer(ndims);
    optimizer.tolerance(1e-6);
    optimizer.maximumIterations(10000);
    optimizer.seedLengths(lens);
    vector<double> final = optimizer.optimize(seeds, fitFunction);
    
    for (uint i=0; i<final.size(); ++i)
      cout << boost::format("%12.8g ") % final[i];
    cout << boost::format("\t\t%.6g\n") % optimizer.finalValue();
    seeds = final;
  }
  
}
