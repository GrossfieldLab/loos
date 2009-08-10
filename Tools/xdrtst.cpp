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
    
    uint block_size = sizeof(loos::internal::XDR::block_type);
    uint offset = 9 * block_size;
    (xf.get())->seekg(offset, ios_base::cur);
    uint nbytes;
    xf.read(&nbytes);
    uint nblocks = nbytes / block_size;
    if (nbytes % block_size != 0)
      ++nblocks;
    offset = nblocks * block_size;
    (xf.get())->seekg(offset, ios_base::cur);
  }
}




int main(int argc, char *argv[]) {
  fstream fs("f.xtc", ios_base::in);
  loos::internal::XDR xfile(&fs);

  scanFile(xfile);
  fs.close();

  cout << "--MARKER--MARKER--MARKER--MARKER--\n";

  // Now try it through trajectory interface...
  loos::AtomicGroup model = loos::createSystem("f.gro");

  cout << model;

  cout << "--MARKER--MARKER--MARKER--MARKER--\n";

  loos::pTraj traj = loos::createTrajectory("f.xtc", model);


  cout << "nframes = " << traj->nframes() << endl;
  cout << "natoms = " << traj->natoms() << endl;

  int n = 0;
  while (traj->readFrame()) {
    cout << format("Frame = %d\n") % n++;
    cout << format("Box = %s\n") % traj->periodicBox();
    traj->updateGroupCoords(model);
    for (uint i=0; i<5; ++i)
      cout << *(model[i]) << endl;
  }

  
}
