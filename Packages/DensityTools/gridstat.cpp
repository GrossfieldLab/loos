/*
  gridstat.cpp

  (c) 2008 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry

      Simple grid statistics...
*/



#include <loos.hpp>
#include <boost/format.hpp>
#include <DensityGrid.hpp>

using namespace std;
using namespace loos;
using namespace loos::DensityTools;




double avgDens(DensityGrid<double>& grid) {
  long i, n = grid.maxGridIndex();
  double avg = 0.0;

  for (i=0; i<n; i++)
    avg += grid(i);

  return(avg/n);
}


double zavgDens(DensityGrid<double>& grid) {
  long n = grid.maxGridIndex();
  long m = 0;
  double avg = 0.0;

  for (long i=0; i<n; ++i) {
    double d = grid(i);
    if (d > 0.0) {
      avg += d;
      ++m;
    }
  }

  avg /= m;
  return(avg);
}


double stdDens(DensityGrid<double>& grid, const double avg) {
  long i, n = grid.maxGridIndex();
  double std = 0.0;

  for (i=0; i<n; i++)
    std += (grid(i) - avg) * (grid(i) - avg);

  return(sqrt(std/(n-1.0)));
}


double zstdDens(DensityGrid<double>& grid, const double avg) {
  long i, n = grid.maxGridIndex();
  double std = 0.0;

  long m = 0;
  for (i=0; i<n; i++)
    if (grid(i) > 0.0) {
      std += (grid(i) - avg) * (grid(i) - avg);
      ++m;
    }

  return(sqrt(std/(m-1.0)));
}


double maxDens(DensityGrid<double>& grid) {
  long i, n = grid.maxGridIndex();
  double max = 0.0;

  for (i=0; i<n; i++)
    if (grid(i) > max)
      max = grid(i);

  return(max);
}



void quickHist(DensityGrid<double>& grid, const double x, const int nbins) {
  long *bins = new long[nbins];
  double delta = x / nbins;
  
  for (int i=0; i<nbins; i++)
    bins[i] = 0;

  long n = grid.maxGridIndex();
  for (long i=0; i<n; i++) {
    int k = static_cast<int>(grid(i) / delta);
    assert(k <= nbins && k >= 0);
    if (k == nbins)
      k = nbins - 1;
    ++bins[k];
  }

  cout << "Quick histogram\n";
  cout << "---------------\n";
  for (int i=0; i<nbins; i++) {
    cout << setprecision(6) << setw(10) << i*delta << "\t" << setprecision(4) << static_cast<double>(i)/nbins << "\t";
    cout << setw(10) << bins[i] << "\t" << setprecision(4) << static_cast<double>(bins[i]) / n << endl;
  }

  delete[] bins;
}


void zAverage(DensityGrid<double>& grid, const int nbins) {
  DensityGridpoint dims = grid.gridDims();

  int chunk_size = dims[2] / nbins;
  long volume = chunk_size * dims[1] * dims[0];

  cout << endl;
  cout << "Z-slice averages\n";
  cout << "----------------\n";
  
  int kk = 0;
  for (int k = 0; k<nbins; k++) {
    // Calculate z-range...
    DensityGridpoint bottom(0,0,k*chunk_size);
    DensityGridpoint top(0,0,chunk_size*(k+1));
    
    GCoord wbottom = grid.gridToWorld(bottom);
    GCoord wtop = grid.gridToWorld(top);

    double avg = 0.0;

    for (int sk = 0; sk < chunk_size && sk+kk < dims[2]; sk++, kk++)
      for (int j=0; j<dims[1]; j++)
	for (int i=0; i<dims[0]; i++)
	  avg += grid(kk, j, i);

    avg /= volume;
    cout << kk << "\t" << wbottom.z() << "\t" << wtop.z() << "\t" << avg << endl;
  }

  if (kk < dims[2]) {
    DensityGridpoint bottom(0,0,kk);
    GCoord wbottom = grid.gridToWorld(bottom);
    double avg = 0.0;

    volume = 0;
    for (; kk < dims[2]; kk++)
      for (int j=0; j<dims[1]; j++)
	for (int i=0; i<dims[0]; i++, volume++)
	  avg += grid(kk, j, i);

    DensityGridpoint top(0,0,kk);
    GCoord wtop = grid.gridToWorld(top);

    avg /= volume;
    cout << kk << "\t" << wbottom.z() << "\t" << wtop.z() << "\t" << avg << endl;
    cout << "Warning- last row adjusted\n";
  }
}




int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "Usage- gridstat bins zbins <file.grid\n";
    exit(-1);
  }

  double nbins = strtod(argv[1], 0);
  double zbins = strtod(argv[2], 0);

  DensityGrid<double> grid;
  cin >> grid;

  cout << "Read in grid of size " << grid.gridDims() << endl;
  cout << "Range is " << grid.minCoord() << " to " << grid.maxCoord() << endl;
  
  double gavg = avgDens(grid);
  double gzavg = zavgDens(grid);
  double gstd = stdDens(grid, gavg);
  double gzstd = zstdDens(grid, gzavg);
  double gmax = maxDens(grid);

  cout << "\n\n* Grid Density Statistics *\n";
  cout << "Grid density is " << gavg << " (" << gstd << ")\n";
  cout << "Grid non-zero avg is " << gzavg << " (" << gzstd << ")\n";
  cout << "Max density is " << gmax << endl << endl;
  quickHist(grid, gmax, nbins);
  zAverage(grid, zbins);


}
