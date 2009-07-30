#include <loos.hpp>
#include <xdr.hpp>

#include <boost/format.hpp>


using namespace std;
using namespace boost;

int main(int argc, char *argv[]) {
  fstream fs("f.xtc", ios_base::in);
  loos::internal::XDR xfile(&fs);

  int natoms;
  int step;
  float time;
  float box[9];

  int magic = 99;
  int r = xfile.read(magic);
  cout << format("magic=%d\n") % magic;
  cout << format("r=%d\n") % r;

  xfile.read(natoms);
  xfile.read(step);
  xfile.read(time);
  xfile.read(box, 9);


  cout << format("natoms=%d, step=%d, time=%f\n") % natoms % step % time;
  cout << "box=(";
  for (int i=0; i<9; ++i)
    cout << box[i] << ((i == 8) ? "" : ",");
  cout << endl;

  vector<float> coords;
  float precision;
  r = loos::internal::xtc::xdrfile_read_compr_coord(coords, &precision, &xfile);
  cout << format("Read of coords: r=%d, prec=%f\n") % r % precision;

  coords.clear();

  r = xfile.read(magic);
  cout << format("magic=%d\n") % magic;
  cout << format("r=%d\n") % r;
  xfile.read(natoms);
  xfile.read(step);
  xfile.read(time);
  xfile.read(box, 9);


  cout << format("natoms=%d, step=%d, time=%f\n") % natoms % step % time;
  cout << "box=(";
  for (int i=0; i<9; ++i)
    cout << box[i] << ((i == 8) ? "" : ",");
  cout << endl;

  r = loos::internal::xtc::xdrfile_read_compr_coord(coords, &precision, &xfile);
  cout << format("Read of coords: r=%d, prec=%f\n") % r % precision;

}
