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


#include <ProgressCounters.hpp>
#include <utils.hpp>
#include <cmath>


namespace loos {

  void SimpleCounter::attach(ObsT* obs) { observers.push_back(obs); }
  void SimpleCounter::detach(ObsT* obs) {
    Observers::iterator i = std::find(observers.begin(), observers.end(), obs);
    if (i == observers.end())
      throw(std::logic_error("Attempting to detach an observer that was never attached"));
    observers.erase(i);
  }

  void SimpleCounter::notify(void) {
    Observers::iterator i;
    for (i = observers.begin(); i != observers.end(); ++i)
      (*i)->update(this);
  }

  void SimpleCounter::finish(void) {
    timer_.stop();
    Observers::iterator i;
    for (i = observers.begin(); i != observers.end(); ++i)
      (*i)->finish(this);
  }

  void SimpleCounter::start(void) {
    count_ = 0;
    timer_.start();
    Observers::iterator i;
    for (i = observers.begin(); i != observers.end(); ++i)
      (*i)->start(this);
  }

  
  uint SimpleCounter::count(void) const { return(count_); }
  double SimpleCounter::elapsed(void) { return(timer_.elapsed()); }

  uint SimpleCounter::remaining(void) { throw(std::logic_error("remaining() is unimplemented")); }

  double SimpleCounter::timeRemaining(void) { throw(std::logic_error("timeRemaining() is unimplemented")); }

  double SimpleCounter::fractionComplete(void) { throw(std::logic_error("fractionComplete() is unimplemented")); }




  // ----------------------------------------------------------

  void EstimatingCounter::setExpected(uint n) { expected = n; }
  uint EstimatingCounter::remaining(void) { return(expected - count_); }
  double EstimatingCounter::timeRemaining(void) {
    double avg = timer_.elapsed() / count_;
    return(remaining() * avg);
  }

  double EstimatingCounter::fractionComplete(void) { return(static_cast<double>(count_) / expected); }


  // ----------------------------------------------------------

  void BasicProgress::start(SimpleCounter*) { os_ << prefix_; }
  void BasicProgress::update(SimpleCounter*) { os_ << msg_; }
  void BasicProgress::finish(SimpleCounter*) { os_ << suffix_; }

  // ----------------------------------------------------------

  void PercentProgress::update(SimpleCounter* s) {
    uint i = static_cast<uint>(floor(s->fractionComplete() * 100.0));
    os_ << i << "% " << msg_ << " (" << timeAsString(s->timeRemaining()) << " remaining)\n";
  }

  void PercentProgress::finish(SimpleCounter* s) {
    os_ << suffix_;
    os_ << "Total elapsed time was " << timeAsString(s->elapsed()) << std::endl;
  }

}

