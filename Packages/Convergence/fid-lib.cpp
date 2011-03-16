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




vecUint findFreeFrames(const vecInt& map) {
  vecUint indices;

  for (uint i=0; i<map.size(); ++i)
    if (map[i] < 0)
      indices.push_back(i);

  return(indices);
}



vecUint assignStructures(AtomicGroup& model, pTraj& traj, const vecUint& frames, const vecGroup& refs) {
  vecUint assignments(frames.size(), 0);
  uint j = 0;

  for (vecUint::const_iterator frame = frames.begin(); frame != frames.end(); ++frame) {
    traj->readFrame(*frame);
    traj->updateGroupCoords(model);

    double mind = numeric_limits<double>::max();
    uint mini = refs.size() + 1;
    
    for (uint i = 0; i < refs.size(); ++i) {
      model.alignOnto(refs[i]);
      double d = model.rmsd(refs[i]);
      if (d < mind) {
        mind = d;
        mini = i;
      }
    }

    assignments[j++] = mini;
  }

  return(assignments);
}


vecUint trimFrames(const vecUint& frames, const double frac) {
  uint bin_size = frac * frames.size();
  uint remainder = frames.size() - static_cast<uint>(bin_size / frac);

  vecUint truncated;
  vecUint::const_iterator vi = frames.begin();
  for (uint i = 0; i<frames.size() - remainder; ++i)
    truncated.push_back(*vi++);
  
  return(truncated);
}




boost::tuple<vecGroup, vecUint> pickFiducials(AtomicGroup& model, pTraj& traj, const vecUint& frames, const double f) {

  // Size of bin
  uint bin_size = f * frames.size();

  // Initialize a RNG to use a uniform random distribution.
  // Use the LOOS generator singleton so the random number stream can
  // be seeded (or automatically seeded) by LOOS...
  boost::uniform_real<> rmap;
  boost::variate_generator< base_generator_type&, boost::uniform_real<> > rng(rng_singleton(), rmap);

  // Vector of AtomicGroup's representing the fiducial structures
  vecGroup fiducials;

  // Track which trajectory frame has been assigned to which fiducial
  vecInt assignments(frames.size(), -1);

  // The indices (frame #'s) of the structures picked to be fiducials
  vecUint refs;
  vecDouble radii;
  
  // Unassigned frames...bootstrap the loop
  vecUint possible_frames = findFreeFrames(assignments);

  // Are there any unassigned frames left?
  while (! possible_frames.empty()) {
    // Randomly pick one
    uint pick = possible_frames[static_cast<uint>(floor(possible_frames.size() * rng()))];

    traj->readFrame(frames[pick]);
    traj->updateGroupCoords(model);

    // Make a copy and assign a new bin # to the fiducial
    AtomicGroup fiducial = model.copy();
    fiducial.centerAtOrigin();
    uint myid = fiducials.size();
    if (assignments[pick] >= 0) {
      cerr << "INTERNAL ERROR - " << pick << " pick was already assigned to " << assignments[pick] << endl;
      exit(-99);
    }

    fiducials.push_back(fiducial);
    refs.push_back(pick);
    
    // Now find the distance from every unassigned frame to this new fiducial (aligning
    // them first), and then sort by distance...
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

    // Pick the first bin_size of them (or however many are remaining)
    // and assign these to the newly picked fiducial
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

  boost::tuple<vecGroup, vecUint> result(fiducials, refs);
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
  vecUint histogram(max+1, 0);

  for (vecInt::const_iterator i = assignments.begin(); i != assignments.end(); ++i)
    if (*i >= 0)
      histogram[*i] += 1;

  return(histogram);
}
