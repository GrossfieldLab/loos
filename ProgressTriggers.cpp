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


#include <loos/ProgressTriggers.hpp>


namespace loos {

  bool TriggerEvery::operator()(SimpleCounter* subj) {
    return(subj->count() % freq == 0);
  }

  void TriggerEvery::setFrequency(const uint i) { freq = i; }


  //----------------------------------------------------------
  void PercentTrigger::setFraction(double frac) { frac_ = frac; }
  void PercentTrigger::reset(void) { chunk_ = 0; }

  bool PercentTrigger::operator()(SimpleCounter* subj) {
    int chunk = static_cast<int>( subj->fractionComplete() / frac_ );
    if (chunk != chunk_) {
      chunk_ = chunk;
      return(true);
    }
    return(false);
  };
  
}
