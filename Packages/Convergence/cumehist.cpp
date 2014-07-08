// Histogram of a time series using an increasingly larger window


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2014, Tod D. Romo and Alan Grossfield
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

using namespace std;
using namespace loos;



vector<double> histogram(const vector<double>& data,
			 const uint nelems,
			 const uint nbins,
			 const double minval,
			 const double maxval) {

  vector<uint> hist(nbins, 0);
  double delta = nbins / (maxval - minval);

  for (uint i=0; i<nelems; ++i) {
    uint bin = (data[i] - minval) * delta;
    if (bin < nbins)
      hist[bin] += 1;
  }
  
  vector<double> h(nbins);
  for (uint i=0; i<nbins; ++i)
    h[i] = static_cast<double>(hist[i]) / nelems;

  return(h);
}


pair<double, double> findMinMax(const vector<double>& data) {
  double max = numeric_limits<double>::min();
  double min = numeric_limits<double>::max();

  for (vector<double>::const_iterator ci = data.begin(); ci != data.end(); ++ci) {
    if (*ci < min)
      min = *ci;
    if (*ci > max)
      max = *ci;
  }

  return(pair<double, double>(min, max));
}


vector<double> readData(const string& fname, const uint col) 
{
  vector< vector<double> > table = readTable<double>(fname);
  vector<double> d(table.size());
  
  for (uint i=0; i<table.size(); ++i)
    d[i] = table[i][col];
  
  return(d);
}



int main(int argc, char *argv[]) {

  if (!(argc == 5 || argc == 7)) {
    cerr << "Usage- " << argv[0] << "datafile col nbins stride [min max]\n";
    exit(-1);
  }

  int k = 1;
  string hdr = invocationHeader(argc, argv);
  string fname(argv[k++]);
  uint col = strtoul(argv[k++], 0, 10);
  uint nbins = strtoul(argv[k++], 0, 10);
  uint stride = strtoul(argv[k++], 0, 10);


  vector<double> data = readData(fname, col);

  double minval, maxval;
  if (k < argc) {
    minval = strtod(argv[k++], 0);
    maxval = strtod(argv[k++], 0);
  } else {
    pair<double,double> r = findMinMax(data);
    minval = r.first;
    maxval = r.second;
    ostringstream oss;
    
    oss << hdr << "\n# min = " << minval << endl << "# max = " << maxval;
    hdr = oss.str();
    
  }

  cout << "# " << hdr << endl;
  
  double factor = (maxval - minval) / nbins;

  for (uint y = stride; y<data.size(); y += stride) {
    vector<double> h = histogram(data, y, nbins, minval, maxval);
    for (uint n=0; n<nbins; ++n) {
      double x = (n + 0.5) * factor + minval;
      cout << x << '\t' << y << '\t' << h[n] << endl;
    }
    cout << endl;
  }
  
}
