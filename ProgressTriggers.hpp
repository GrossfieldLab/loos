/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008-2009, Tod D. Romo, Alan Grossfield
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

#if !defined(PROGRESSTRIGGERS_HPP)
#define PROGRESSTRIGGERS_HPP

#include <loos_defs.hpp>
#include <ProgressCounters.hpp>


namespace loos {


  class TriggerEvery {
  public:
    TriggerEvery(const uint i) : freq(i) { }
    bool operator()(SimpleCounter* subj) {
      return(subj->count() % freq == 0);
    }

    void setFrequency(const uint i) { freq = i; }

  private:
    uint freq;
  };



  class PercentTrigger {
  public:
    PercentTrigger(double frac) : frac_(frac), chunk_(0) { }
    void setFraction(double frac) { frac_ = frac; }
    void reset(void) { chunk_ = 0; }

    bool operator()(SimpleCounter* subj) {
      int chunk = static_cast<int>( subj->fractionComplete() / frac_ );
      if (chunk != chunk_) {
        chunk_ = chunk;
        return(true);
      }
      return(false);
    };

  private:
    double frac_;
    int chunk_;
  };


}


#endif
