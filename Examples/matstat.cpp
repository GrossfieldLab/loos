#include <loos.hpp>
#include <MatrixReader.hpp>

typedef unsigned int uint;


int main(int argc, char *argv[]) {
  RawAsciiReader<float> reader;
  
  RawAsciiReader<float>::Result res = reader.read(argv[1]);
  uint m = boost::get<1>(res);
  uint n = boost::get<2>(res);
  float *p = boost::get<0>(res);

  cout << "Read in a " << m << "x" << n << " matrix from " << argv[1] << endl;
  float mean = 0.0, max=p[0], min=p[0];
  uint k, kk = m*n;

  for (k=0; k<kk; k++) {
    mean += p[k];
    if (p[k] > max)
      max = p[k];
    if (p[k] < min)
      min = p[k];
  }

  mean /= kk;

  cout << "Max = " << max << endl;
  cout << "Min = " << min << endl;
  cout << "Avg = " << mean << endl;

}

