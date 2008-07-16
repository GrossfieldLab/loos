#include <iostream>
#include <string>
#include <loos.hpp>


int main(int argc, char *argv[]) {
  Amber file(argv[1]);
  AmberTraj at(argv[2], file.size());
  DCDWriter dcd(argv[3]);
  dcd.setHeader(file.size(), at.nframes(), 1e-3, at.hasPeriodicBox());
  dcd.writeHeader();

  while (at.readFrame()) {
    at.updateGroupCoords(file);
    dcd.writeFrame(file);
  }
}

