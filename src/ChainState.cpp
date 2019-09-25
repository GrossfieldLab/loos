

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2019, Alan Grossfield
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

#include "ChainState.hpp"

namespace loos {

    void ChainState::computeChainState(const AtomicGroup &group,
                                       const GCoord &normal,
                                       StateVector &segs) {
        for (uint i=0; i< _num_segs; ++i) {
            loos::GCoord vec = group[i]->coords() - group[i+1]->coords();
            double cosine = vec*normal/vec.length();
            cosine = fmin(cosine, 1.0);
            cosine = fmax(cosine, -1.0);
            //std::cerr << i << "  " << cosine << "  " << _bin_width << "  ";
            //std::cerr << vec << "  " << normal << "  ";
            segs[i] = static_cast<int>((cosine) / _bin_width);
            //std::cerr << segs[i] << std::endl;
        }

        // Store the state of this chain
        state_counts[segs] += 1;
        counts++;
    }

    void ChainState::computeChainState(const AtomicGroup &group,
                                      const GCoord &normal) {
    StateVector segs = StateVector(_num_segs);
    computeChainState(group, normal, segs);
    }

    double ChainState::getStateProb(const StateVector &segs) {
        if (state_counts.count(segs)) {
            return static_cast<double>(state_counts[segs])/counts;
        }
        else {
            return 0.0;
        }
    }

    std::set<std::pair<StateVector, uint>, Comparator > ChainState::getAllProbs() {
        // Copy the elements of the map into a set
        std::set<std::pair<StateVector, uint>, Comparator >
            state_set(state_counts.begin(),
                      state_counts.end(),
                      compareStateProbs);
        return(state_set);
    }

    double ChainState::entropy() {
        double ent = 0.0;
        for (std::map<StateVector, uint>::iterator s = state_counts.begin();
                                                  s!= state_counts.end();
                                                  ++s) {
            double prob = static_cast<double>(s->second) / num_counts();
            ent -= prob * log(prob);
            }
        return ent;
    }

    double ChainState::relative_entropy(const std::map<StateVector, double> &ref) {
        /**
        The map value is double, not uint, because it presumably will have been
        read in from a file that was already normalized.
        */
        double ent = 0.0;
        for (std::map<StateVector, uint>::iterator s = state_counts.begin();
                                                   s!= state_counts.end();
                                                   ++s) {
            if (ref.count(s->first)) {
                double p = static_cast<double>(s->second) / num_counts();
                double ratio = p / ref.at(s->first);
                ent += p * log(ratio);
                }
            }
        return(ent);
    }

    RefChainDist::RefChainDist(const std::string &filename) {
        readInput(filename);
    }

    void RefChainDist::readInput(const std::string &filename) {
        std::ifstream ifs(filename.c_str());
        if (!ifs) {
            throw(FileOpenError(filename, "Couldn't open reference distribution file"));
        }
        std::string input;
        bool first = true;
        int prev_size = 0;
        while (std::getline(ifs, input)) {
            // Skip blank lines and lines starting with "#"
            if ( (input.length() == 0) || (input[0] == '#' ) ) {
                // do nothing
                continue;
            }
            std::istringstream ist(input);

            StateVector state_vector;
            double prob;
            int state;
            ist >> prob;
            while (ist.good()) {
                ist >> state;
                state_vector.push_back(state);
            }
            if (first) {
                first = false;
                prev_size = state_vector.size();
            }
            else {
                if (state_vector.size() != prev_size) {
                    throw(LOOSError("all reference states must be same length"));
                }
            }
            state_dist[state_vector] = prob;
        }

    }

    double RefChainDist::relative_entropy(RefChainDist &ref) {

        double ent = 0.0;
        for (std::map<StateVector, double>::iterator s = state_dist.begin();
                                                     s!= state_dist.end();
                                                     ++s) {
            if (ref.state_dist.count(s->first)) {
                double p = s->second;
                double ratio = p / ref.state_dist.at(s->first);
                ent += p * log(ratio);
                }
            }
        return(ent);
    }


}
