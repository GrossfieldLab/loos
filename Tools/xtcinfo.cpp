#include <loos.hpp>

using namespace std;
using namespace loos;


const int magic = 1995;


struct Header {
  uint natoms, step;
  float time, box[9];
};



bool readFrameHeader(internal::XDRReader& xdr, Header& hdr) {
  int magic_no;
  int ok = xdr.read(magic_no);
  if (!ok)
    return(false);
  if (magic_no != magic) {
    cerr << "Invalid XTC magic number (got " << magic_no << " but expected " << magic << ")";
    exit(-1);
  }

  // Defer error-checks until the end...
  xdr.read(hdr.natoms);

  xdr.read(hdr.step);
  xdr.read(hdr.time);
  ok = xdr.read(hdr.box, 9);
  if (!ok)
    return(false);

  uint block_size = sizeof(internal::XDRReader::block_type);

  size_t offset = 0;
  uint nbytes = 0;

  if (hdr.natoms <= 9) {
    nbytes = hdr.natoms * 3 * sizeof(float);
    uint dummy;
    xdr.read(dummy);
  } else {
    offset = 9 * block_size;
    xdr.get()->seekg(offset, ios_base::cur);
    xdr.read(nbytes);
  }

  uint nblocks = nbytes / block_size;
  if (nbytes % block_size != 0)
    ++nblocks;

  offset = nblocks * block_size;
  xdr.get()->seekg(offset, ios_base::cur);

  return(true);
}





int main(int argc, char* argv[]) {
  if (argc != 2) {
    cerr << "Usage- xtcinfo filename\n";
    exit(-1);
  }

  ifstream ifs(argv[1]);
  internal::XDRReader xdr(&ifs);
  uint frameno = 0;
  Header hdr;

  cout << boost::format("#%8s %10s %10s %10s\n")
    % "Frame"
    % "NAtoms"
    % "Step"
    % "Time";


  while (readFrameHeader(xdr, hdr)) {
    cout << boost::format(" %8d %10d %10d %10.1f\n")
      % frameno
      % hdr.natoms
      % hdr.step
      % hdr.time;

    ++frameno;
  }


}
