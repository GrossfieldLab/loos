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


namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

// @cond TOOLS_INTERNAL 


string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\tCumulative or windowed histogram\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tThis tool can calculate either a cumulative or a windowed histogram.  The former is\n"
    "made by calculating the histogram of the input data up to time t.  This is written out\n"
    "as a row of data suitable for plotting with gnuplot using the splot command.  Each row\n"
    "then corresponds to calculating the histogram with more points.  The alternative, is\n"
    "to only calculate the histogram over a window that is slid along the data.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\tchist torsion_data >torsion_hist.asc\n"
    "This example uses the defaults, which assumes the column to histogram is column 1\n"
    "(i.e. the second column, since colunm indices are 0-based), wich 20 bins, a stride\n"
    "through the data of 10 (every 10th datapoint is used in the histogram), the range\n"
    "of the histogram is automatically determined from the data, and the histogram type\n"
    "is cumulative.\n"
    "\n"
    "\tchist --min -180 --max 180 --nbins 50 --stride 2 torsion_data >torsion_hist.asc\n"
    "This example is similar to the previous, except that the histogram range is explicitly\n"
    "set to -180 to 180, 50 bins are used, and every other datapoint is taken.\n"
    "\n"
    "\tchist --mode window --window 250 torsion_data.asc >torsion_hist.asc\n"
    "This example calculates a windowed histogram using 250 datapoints per histogram,\n"
    "each window is slid 10 points down (the default for --stride)\n\n";

  return(msg);
}


struct ToolOptions : public opts::OptionsPackage {
public:
  enum ToolMode { CUMULATIVE, WINDOW };
  static const unsigned char SETMIN_F = 0x01;
  static const unsigned char SETMAX_F = 0x02;

  ToolOptions() : col(1),
                  nbins(20),
                  window(100),
                  stride(10),
                  mode_string("cume"),
                  minval(0),
                  maxval(0),
                  range_flags(0),
                  mode(CUMULATIVE)
    {}

  


  void addGeneric(po::options_description& o) {
    o.add_options()
      ("column,C", po::value<uint>(&col)->default_value(col), "Data column to use")
      ("nbins,N", po::value<uint>(&nbins)->default_value(nbins), "Number of bins in histogram")
      ("window", po::value<uint>(&window)->default_value(window), "Histogram window size")
      ("stride", po::value<uint>(&stride)->default_value(stride), "Stride through trajectory for cumulative histogram mode, or how far to slide the window")
      ("mode", po::value<string>(&mode_string)->default_value(mode_string), "Histogram mode: cume or window")
      ("min", po::value<double>(), "Set min value for histogram range")
      ("max", po::value<double>(), "Set max value for histogram range");
  }

  bool postConditions(po::variables_map& map) {
    if (mode_string == "cume")
      mode = CUMULATIVE;
    else if (mode_string == "window")
      mode = WINDOW;
    else {
      cerr << "ERROR- '" << mode_string << "' is an unknown mode.  Must be either 'cume' or 'window'\n";
      return(false);
    }

    if (map.count("min")) {
      minval = map["min"].as<double>();
      range_flags |= SETMIN_F;
    }
    if (map.count("max")) {
      maxval = map["max"].as<double>();
      range_flags |= SETMAX_F;
    }
    return(true);
  }


  string print() const {
    ostringstream oss;

    oss << boost::format("col=%d,nbins=%d,window=%d,stride=%d,mode='%s'")
      % col
      % nbins
      % window
      % stride
      % mode_string;
    return(oss.str());
  }
  
  

  uint col;
  uint nbins;
  uint window;
  uint stride;
  string mode_string;
  double minval, maxval;
  unsigned char range_flags;
  ToolMode mode;
};

// @endcond


vector<double> histogram(const vector<double>& data,
                         const uint start,
			 const uint end,
			 const uint nbins,
			 const double minval,
			 const double maxval) {

  vector<uint> hist(nbins, 0);
  double delta = nbins / (maxval - minval);
  
  for (uint i=start; i<end; ++i) {
    uint bin = (data[i] - minval) * delta;
    if (bin < nbins)
      hist[bin] += 1;
  }
  
  vector<double> h(nbins);
  uint nelems = end - start;
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

  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  ToolOptions* topts = new ToolOptions;
  opts::RequiredArguments* ropts = new opts::RequiredArguments("datafile", "Name of file to histogram");

  opts::AggregateOptions options;
  options.add(bopts).add(topts).add(ropts);
  if (!options.parse(argc, argv))
    exit(-1);


  vector<double> data = readData(ropts->value("datafile"), topts->col);

  double minval, maxval;
  pair<double,double> r = findMinMax(data);
  minval = r.first;
  maxval = r.second;
  
  if (topts->range_flags & ToolOptions::SETMIN_F)
    minval = topts->minval;
  if (topts->range_flags & ToolOptions::SETMAX_F)
    maxval = topts->maxval;

  cout << "# " << hdr << endl;
  cout << "# min = " << minval << endl << "# max = " << maxval << endl;
  
  double factor = (maxval - minval) / topts->nbins;

  if (topts->mode == ToolOptions::CUMULATIVE) {

    for (uint y = topts->stride; y<data.size(); y += topts->stride) {
      vector<double> h = histogram(data, 0, y, topts->nbins, minval, maxval);
      for (uint n=0; n<topts->nbins; ++n) {
        double x = (n + 0.5) * factor + minval;
        cout << x << '\t' << y << '\t' << h[n] << endl;
      }
      cout << endl;
    }

  } else {

    for (uint y = 0; y<data.size() + topts->window; y += topts->stride) {
      vector<double> h = histogram(data, y, y+topts->window, topts->nbins, minval, maxval);
      for (uint n=0; n<topts->nbins; ++n) {
        double x = (n + 0.5) * factor + minval;
        cout << x << '\t' << y << '\t' << h[n] << endl;
      }
      cout << endl;
      
    }
    
  }
    
}
