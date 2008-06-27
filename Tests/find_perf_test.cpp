/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo
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


#define ERROR_CHECK


int main(int argc, char *argv[]) {
  long n = atol(argv[1]);
  PDB pdb(argv[2]);

  int a = pdb.minId();
  int b = pdb.maxId();
  int d = b-a;

  int *indices = new int[n];
  int i;

  cerr << "Generating indices...\n";
  srand48(time(0));
  for (i=0; i<n; i++)
    indices[i] = (int)(drand48() * d) + a;

  cerr << "Searching...";
  pAtom pa;

  for (i=0; i<n; i++) {
    pa = pdb.findById(indices[i]);
#if defined(ERROR_CHECK)
    if (pa->id() != indices[i]) {
      cerr << "***ERROR*** at index " << i << "(" << indices[i] << ")\n";
      exit(-20);
    } else if (pa == 0)
      exit(-10);
#endif
  }

  cerr << "done\n";
}
