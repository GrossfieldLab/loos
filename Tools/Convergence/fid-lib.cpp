/*
  fid-lib
  
  Fiducial library
*/


/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2010, Tod D. Romo
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



#include "fid-lib.hpp"

using namespace std;
using namespace loos;




vecUint findFreeFrames(const vector<int>& map) {
  vecUint indices;

  for (uint i=0; i<map.size(); ++i)
    if (map[i] < 0)
      indices.push_back(i);

  return(indices);
}






boost::tuple<vecInt, vecUint, vecGroup, vecDouble> assignFrames(AtomicGroup& model, pTraj& traj, const vecUint& frames, const double f) {

  uint bin_size = f * frames.size();
  boost::uniform_real<> rmap;
  boost::variate_generator< base_generator_type&, boost::uniform_real<> > rng(rng_singleton(), rmap);

  vecGroup fiducials;
  vecInt assignments(frames.size(), -1);
  vecUint refs;
  vecDouble radii;
  
  vecUint possible_frames = findFreeFrames(assignments);
  while (! possible_frames.empty()) {
    uint pick = possible_frames[static_cast<uint>(floor(possible_frames.size() * rng()))];

    traj->readFrame(frames[pick]);
    traj->updateGroupCoords(model);

    AtomicGroup fiducial = model.copy();
    fiducial.centerAtOrigin();
    uint myid = fiducials.size();
    if (assignments[pick] >= 0) {
      cerr << "INTERNAL ERROR - " << pick << " pick was already assigned to " << assignments[pick] << endl;
      exit(-99);
    }

    fiducials.push_back(fiducial);
    refs.push_back(pick);
    
    vector<double> distances(assignments.size(), numeric_limits<double>::max());
    for (uint i = 0; i<assignments.size(); ++i) {
      if (assignments[i] >= 0)
        continue;
      traj->readFrame(frames[i]);
      traj->updateGroupCoords(model);
      model.centerAtOrigin();
      model.alignOnto(fiducial);
      distances[i] = model.rmsd(fiducial);
    }

    vecUint indices = sortedIndex(distances);
    uint picked = 0;
    double maxd = 0.0;
    for (uint i=0; i<assignments.size() && picked < bin_size; ++i) {
      if (assignments[indices[i]] < 0) {
        assignments[indices[i]] = myid;
        ++picked;
        if (distances[indices[i]] > maxd)
          maxd = distances[indices[i]];
      }
    }
    radii.push_back(maxd);
    possible_frames = findFreeFrames(assignments);
  }

  // Safety check...
  for (vecInt::const_iterator i = assignments.begin(); i != assignments.end(); ++i)
    if (*i < 0)
      throw(runtime_error("A frame was not assigned in binFrames()"));

  boost::tuple<vecInt, vecUint, vecGroup, vecDouble> result(assignments, refs, fiducials, radii);
  return(result);
}




int findMaxBin(const vecInt& assignments) {
  
  int max_bin = -1;
  for (vecInt::const_iterator i = assignments.begin(); i != assignments.end(); ++i)
    if (*i > max_bin)
      max_bin = *i;

  return(max_bin);
}



vecUint histogramBins(const vecInt& assignments) {
  int max = findMaxBin(assignments);
  vecUint histogram(max, 0);

  for (vecInt::const_iterator i = assignments.begin(); i != assignments.end(); ++i)
    if (*i >= 0)
      histogram[*i] += 1;

  return(histogram);
}
