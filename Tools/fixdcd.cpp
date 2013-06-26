
/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2013, Tod D. Romo, Grossfield Lab
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester

  This package (LOOS) is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation under version 3 of the License.

  This package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/




#include <loos.hpp>

using namespace loos;
using namespace std;






int scanTrajectory(const char* fname, bool* swabbing)
{

    DCD dcd(fname);
    *swabbing = !dcd.nativeFormat();
    

    int n = 0;
    while (dcd.readFrame())
        ++n;

    if (n != static_cast<int>(dcd.nframes())) {
        cout << boost::format("%s claims to have %d frames.\n") % fname % dcd.nframes();
        cout << boost::format("--> Scanning found %d frames.\n") % n;
    } else
        n = -1;   // Signal no-update
    
    return(n);
}



void fixDCD(const char* fname) 
{
    bool swabbing;

    int nframes = scanTrajectory(fname, &swabbing);
    if (nframes < 0)
        return;
    

    
    fstream file(fname, ios_base::in | ios_base::out);

    if (file.fail()) {
        cerr << "Error- cannot open " << fname << endl;
        exit(-2);
    }

    if (swabbing)
        nframes = swab(nframes);
    
    file.seekp(8);
    file.write(reinterpret_cast<char*>(&nframes), sizeof(nframes));

    file.seekp(20);
    file.write(reinterpret_cast<char*>(&nframes), sizeof(nframes));

}



int main(int argc, char *argv[]) 
{

    if (argc < 2) {
        cerr << "Usage- fixdcd dcdfile [dcdfile ...]\n";
        exit(-1);
    }

    DCD::setSuppression(true);
    
    for (int k=1; k<argc; ++k)
        fixDCD(argv[k]);

}

        
    

