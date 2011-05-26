/*
  Gridgauss

  (c) 2009 Tod D. Romo
*/


#include <loos.hpp>
#include <boost/format.hpp>
#include <sgrid.hpp>

using namespace std;
using namespace loos;
using namespace lab;


SGrid<double> buildGaussian3d(const int w, const double sigma) {
  SGrid<double> kernel(GCoord(0,0,0),GCoord(0,0,0),w);

  double c = static_cast<double>(w-1)/2.0;
  double a = 1.0 / (sigma * sqrt(2*PI));
  double b = - 1.0 / (2*sigma*sigma);

  double sum = 0.0;

  for (int k=0; k<w; ++k) {
    double z = k-c;
    z *= z;

    for (int j=0; j<w; ++j) {
      double y = j-c;
      y *= y;
      
      for (int i=0; i<w; ++i) {
        double x = i-c;
        x *= x;
        double r = sqrt(x + y + z);
        double val = a * exp(b * r * r);

        kernel(k,j,i) = val;
        sum += val;
      }
    }
  }

  for (long i=0; i<w*w*w; ++i)
    kernel(i) /= sum;

  return(kernel);
}



void showGrid(const SGrid<double>& g, const string& msg) {
  cerr << msg << endl;
  SGridpoint n = g.gridDims();
  for (int k=0; k<n.z(); ++k) {
    for (int j=0; j<n.y(); ++j) {
      for (int i=0; i<n.x(); ++i)
        cerr << boost::format("%12.6e ") % g(k, j, i);
      cerr << endl;
    }
    cerr << endl;
  }

}


int main(int argc, char *argv[]) {

  if (argc != 3) {
    cerr << "Usage- gridgauss width sigma <grid >output\n";
    exit(-1);
  }
  
  string hdr = invocationHeader(argc, argv);
  int width = atoi(argv[1]);
  double sigma = strtod(argv[2], 0);

  //  if (width % 2 == 0) {
  //    cerr << "Error- width must be an odd number.\n";
  //    exit(-1);
  //  }

  SGrid<double> kernel = buildGaussian3d(width, sigma);
  showGrid(kernel, "Kernel");

  SGrid<double> grid;
  cin >> grid;
  SGridpoint dims = grid.gridDims();

  SGrid<double> convolved(grid);
  convolved.zero();

  int kn = dims[2];
  int jn = dims[1];
  int in = dims[0];

  int h = (width - 1) / 2;

  for (int k=h; k<kn-h; ++k)
    for (int j=h; j<jn-h; ++j)
      for (int i=h; i<in-h; ++i) {
        
        double val = 0.0;
        for (int c=0; c<width; ++c)
          for (int b=0; b<width; ++b)
            for (int a=0; a<width; ++a)
              val += grid(k+c-h, j+b-h, i+a-h) * kernel(c, b, a);

        convolved(k, j, i) = val;
      }


  cout << convolved;
}


  
