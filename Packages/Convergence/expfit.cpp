/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2011, Tod D. Romo
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester

  This package (LOOS) is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation under version 3 of the License.

  This package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/



#include <loos.hpp>
#include <Simplex.hpp>
#include <limits>

using namespace std;
using namespace loos;


typedef pair<double, double>     LPoint;
typedef vector<double>           vecDouble;
typedef vector<vecDouble>        vvecDouble;

struct ExponentialFit {
  ExponentialFit(const vector<LPoint>& datapoints) : _datapoints(datapoints) { }

  double operator()(const vector<double>& v) {
    double sum = 0.0;

    for (vector<double>::const_iterator i = v.begin(); i != v.end(); ++i)
      if (*i < 0.0)
        return(numeric_limits<double>::max());

    for (vector<LPoint>::iterator j = _datapoints.begin(); j != _datapoints.end(); ++j) {
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


  vector<LPoint> _datapoints;
};


vvecDouble readData(const string& fname) {
  ifstream ifs(fname.c_str());
  if (!ifs) {
      cerr << "Error- unable to open " << fname << endl;
      exit(-1);
  }
  
  return(readTable<double>(ifs));
}



string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\tExponential fit for bootstrapped-bcom/bcom output\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tThis tool calculates a multi-exponential fit to the bootstrapped-bcom/bcom data.\n"
    "A Nelder-Mead Simplex is used as the optimizer.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\texpfit bcom.asc boot_bcom.asc 5 2 0.7 10 0.3 100\n"
    "This tries to fit bcom.asc and boot_bcom.asc using 5 replicas and using a 2-exponential\n"
    "with inital coefficients of 0.7 and 0.3 and initial correlation times of 10 and 100\n"
    "respectively.\n"
    "\n"
    "SEE ALSO\n"
    "\tbcom, boot_bcom, bootstrap_overlap.pl\n";

  return(msg);
}



int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  if (argc < 5) {
    cerr << "Usage- expfit bcom.asc boot_bcom.asc nreps ndims constant-1 time-1 [constant-2 time-2 ...]\n";
    cerr << fullHelpMessage();
    exit(-1);
  }


  int k = 1;
  vvecDouble bcom = readData(argv[k++]);
  vvecDouble bbcom = readData(argv[k++]);

  if (bcom.size() != bbcom.size()) {
    cerr << boost::format("Error- bcom has %d datapoints but bbcom has %d\n") % bcom.size() % bbcom.size();
    exit(-10);
  }

  vector<LPoint> datapoints;
  for (uint i=0; i<bcom.size(); ++i)
    datapoints.push_back(LPoint(bcom[i][0], bbcom[i][1]/bcom[i][1]));

  int nreps = atoi(argv[k++]);
  int ndims = atoi(argv[k++]);
  vector<double> seeds;
  vector<double> lens;
  ndims *= 2;
  if (argc - 5 != ndims) {
    cerr << boost::format("Error- only %d seeds were specified, but require %d for %d dimensions\n")
      % (argc - 5)
      % (ndims)
      % (ndims/2);
    exit(-1);
  }
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
