#include <loos.hpp>
#include <boost/format.hpp>

using namespace boost;

int main(int argc, char *argv[]) {
  CCPDB ccpdb(argv[1]);
  
  cout << format("There are %u frames with %u atoms.\n") % ccpdb.nframes() % ccpdb.natoms();
  int i = 0;
  while (ccpdb.readFrame()) {
    cout << format("Reading frame %d...\n") % i++;
  }

  cout << format("Read in a total of %d frames.\n") % i;
}
