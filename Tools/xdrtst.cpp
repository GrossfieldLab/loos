#include <loos.hpp>
#include <xdr.hpp>

#include <boost/format.hpp>


using namespace std;
using namespace boost;


bool readHeader(loos::internal::XDR& xf) {
  int natoms;
  int step;
  float time;
  float box[9];

  int magic = 99;
  int r = xf.read(magic);
  if (!r)
    return(false);
  cout << format("magic=%d\n") % magic;

  xf.read(natoms);
  xf.read(step);
  xf.read(time);
  xf.read(box, 9);

  cout << format("natoms=%d, step=%d, time=%f\n") % natoms % step % time;
  cout << "box=(";
  for (int i=0; i<9; ++i)
    cout << box[i] << ((i == 8) ? ")\n" : ",");

  return(true);
}


void scanFile(loos::internal::XDR& xf) {
 
  for (int n = 0; readHeader(xf); ++n) {
    cout << "-- FRAME #" << n << endl;
    
    uint offset = 9 * xf.block_size();
    (xf.get())->seekg(offset, ios_base::cur);
    uint nbytes;
    xf.read(&nbytes);
    uint nblocks = nbytes / xf.block_size();
    if (nbytes % xf.block_size() != 0)
      ++nblocks;
    offset = nblocks * xf.block_size();
    (xf.get())->seekg(offset, ios_base::cur);
  }
}




int main(int argc, char *argv[]) {
  fstream fs("f.xtc", ios_base::in);
  loos::internal::XDR xfile(&fs);

  scanFile(xfile);
  fs.close();

  // Now try it through trajectory interface...
  loos::XTC<float> xtc("f.xtc");
  cout << "nframes = " << xtc.nframes() << endl;
  cout << "natoms = " << xtc.natoms() << endl;

  int n = 0;
  while (xtc.readFrame()) {
    cout << format("Frame = %d\n") % n++;
    cout << format("Box = %s\n") % xtc.periodicBox();
    cout << "First 5 coords:\n";
    vector<loos::GCoord> crds = xtc.coords();
    for (int i=0; i<5; ++i)
      cout << "\t" << crds[i] << endl;
    cout << endl;
  }
  
}
