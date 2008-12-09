#include <iostream>
#include <string>
#include <loos.hpp>

using namespace std;
using namespace loos;


int main(int argc, char *argv[]) {
  Amber file(argv[1]);
  AmberTraj at(argv[2], file.size());
  DCDWriter dcd(argv[3]);
  dcd.setHeader(file.size(), at.nframes(), 1e-3, at.hasPeriodicBox());
  dcd.writeHeader();

  uint i=0;
  cout << "There are " << file.size() << " atoms and " << at.nframes() << " frames.\n";
  while (at.readFrame()) {
    cerr << "Processing frame " << i++ << endl;
    at.updateGroupCoords(file);
    dcd.writeFrame(file);
  }
}

