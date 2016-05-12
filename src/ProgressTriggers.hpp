/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo, Alan Grossfield
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

#if !defined(LOOS_PROGRESSTRIGGERS_HPP)
#define LOOS_PROGRESSTRIGGERS_HPP

#include <loos_defs.hpp>
#include <ProgressCounters.hpp>


namespace loos {

  //! Trigger every i-iterations
  class TriggerEvery {
  public:
    TriggerEvery(const uint i) : freq(i) { }
    bool operator()(SimpleCounter*);
    void setFrequency(const uint);

  private:
    uint freq;
  };


  //! Trigger whenever at least frac percent more iterations have happened
  /**
   * This trigger tracks which fractional update it's in.  If you want
   * to reuse the trigger, you must reset it prior to starting up the
   * counter again...
   */
  class PercentTrigger {
  public:
    PercentTrigger(double frac) : frac_(frac), chunk_(0) { }
    void setFraction(double);
    void reset(void);

    bool operator()(SimpleCounter*);

  private:
    double frac_;
    int chunk_;
  };


}


#endif
