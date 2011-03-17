/*
  gridtst

  (c) 2009 Tod D. Romo
*/



#include <loos.hpp>
#include <sgrid.hpp>


using namespace std;
using namespace loos;
using namespace lab;


int main(int argc, char *argv[]) {
  int n = atoi(argv[1]);

  SGrid<float> grid(GCoord(0,0,0), GCoord(n,n,n), n);

  int x=0;
  for (int k=0; k<n; ++k)
    for (int j=0; j<n; ++j)
      for (int i=0; i<n; ++i)
        grid(k,j,i) = x++;

  cout << grid;
}
