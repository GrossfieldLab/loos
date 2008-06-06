/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
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
#include <AtomicGroup.hpp>
#include <ensembles.hpp>


// Assume all groups are already sorted or matched...

AtomicGroup loos::averageStructure(const vector<AtomicGroup>& ensemble) {
  AtomicGroup avg = ensemble[0].copy();

  // First, zap our coords...
  int n = avg.size();
  int i;
  for (i=0; i<n; i++)
    avg[i]->coords() = GCoord(0.0, 0.0, 0.0);

  // Now, accumulate...
  vector<AtomicGroup>::const_iterator j;
  for (j = ensemble.begin(); j != ensemble.end(); j++) {
    for (i = 0; i<n; i++)
      avg[i]->coords() += (*j)[i]->coords();
  }

  for (i=0; i<n; i++)
    avg[i]->coords() /= ensemble.size();

  return(avg);
}



boost::tuple<vector<XForm>,greal,int> loos::iterativeAlignment(vector<AtomicGroup>& ensemble, greal threshold, int maxiter) {
  int iter = 0;
  int n = ensemble.size();
  greal rms;
  vector<XForm> xforms(n);
  AtomicGroup avg;
  AtomicGroup target = averageStructure(ensemble);

  do {
    for (int i = 0; i<n; i++) {
      GMatrix M = ensemble[i].alignOnto(target);
      xforms[i].concat(M);
    }

    avg = averageStructure(ensemble);
    rms = avg.rmsd(target);
    target = avg;

#if defined(DEBUG)
    cerr << "loos::iterativeAlignment - iter = " << iter << ", rms = " << rms << endl;
#endif

  } while (rms > threshold && ++iter <= maxiter);
  
  boost::tuple<vector<XForm>, greal, int> res(xforms, rms, iter);
  return(res);
}
