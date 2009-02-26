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




#if !defined(PROGRESSCOUNTERS_HPP)
#define PROGRESSCOUNTERS_HPP


#include <iostream>
#include <string>
#include <stdexcept>

#include <loos_defs.hpp>
#include <timer.hpp>


namespace loos {

  class SimpleCounter;

  class AbstractObserver {
  public:
    virtual ~AbstractObserver() { }
    
    virtual void start(SimpleCounter*) =0;
    virtual void finish(SimpleCounter*) =0;
    virtual void update(SimpleCounter*) =0;
  };


  //! Basic progress counter object, defining the interface...
  /**
   * The idea here is that SimpleCounter and its children are
   * "observable" objects (a la the Observer pattern).  The
   * SimpleCounter handles the basic time-keeping functions and
   * forwards messages along to its observers when certain events
   * happen, such as starting, stopping, and updating.  Not all
   * functions are implemented though.  For example, the SimpleTimer
   * doesn't have an estimate for the total number of updates it will
   * receive, so the estimating functions will throw an error.
   *
   * This class is not actually meant to be used by itself, but via
   * the aggregator ProgressCounter class defined later...
   */
  class SimpleCounter  {
    typedef AbstractObserver ObsT;
    typedef std::vector<ObsT*> Observers;
  public:
  
    SimpleCounter() : count_(0) { }
    virtual ~SimpleCounter() { }

    void attach(ObsT* obs) { observers.push_back(obs); }
    void detach(ObsT* obs) {
      Observers::iterator i = std::find(observers.begin(), observers.end(), obs);
      if (i == observers.end())
        throw(std::logic_error("Attempting to detach an observer that was never attached"));
      observers.erase(i);
    }

    //! Notify observers that an update should occur
    virtual void notify(void) {
      Observers::iterator i;
      for (i = observers.begin(); i != observers.end(); ++i)
        (*i)->update(this);
    }

    //! Notify observers that we've finished with our calculation
    virtual void finish(void) {
      timer_.stop();
      Observers::iterator i;
      for (i = observers.begin(); i != observers.end(); ++i)
        (*i)->finish(this);
    }

    //! Notify observers that we're starting a calculation
    virtual void start(void) {
      count_ = 0;
      timer_.start();
      Observers::iterator i;
      for (i = observers.begin(); i != observers.end(); ++i)
        (*i)->start(this);
    }


    //! Number of iterations we've seen so far
    uint count(void) const { return(count_); }

    //! Total elapsed wall-time
    virtual double elapsed(void) { return(timer_.elapsed()); }

    //! Remaining iterations (if applicable)
    virtual uint remaining(void) { throw(std::logic_error("remaining() is unimplemented")); }

    //! Remaining time (estimated, again if applicable)
    virtual double timeRemaining(void) { throw(std::logic_error("timeRemaining() is unimplemented")); }

    //! Percent complete (if applicable)
    virtual double fractionComplete(void) { throw(std::logic_error("fractionComplete() is unimplemented")); }

  protected:
    uint count_;
    Timer<WallTimer> timer_;
    Observers observers;
  };


  //! A progress counter that can estimate how much time is left
  /**
   * This class has to know the count of expected iterations so it can
   * estimate how much time is left.  It is assumed that the execution
   * time per iteration is largely constant.
   */
  class EstimatingCounter : public SimpleCounter {
  public:
    EstimatingCounter(const uint n) : expected(n) { }

    //! Alter the expected count
    void setExpected(uint n) { expected = n; }

    //! Returns the number of iterations left
    uint remaining(void) { return(expected - count_); }

    //! Estimates the amount of time left using the current average
    //! time per iteration
    double timeRemaining(void) {
      double avg = timer_.elapsed() / count_;
      return(remaining() * avg);
    }

    //! Returns the percent completed so far
    double fractionComplete(void) { return(static_cast<double>(count_) / expected); }

  protected:
    uint expected;
  };


  //! This is a simple "trigger" for use as a default
  class TriggerAlways {
  public:
    bool operator()(SimpleCounter*) { return(true); }
  };


  //! The progress counter front-end
  /**
   *
   * This class combines a counter with a trigger.  You can pick
   * whether or not you use a simple counter or an estimating counter,
   * for example, and then combine that with a criterion for firing
   * off a message to any observers that they need to update their
   * display/output.  The trigger is just a functor that returns a
   * true value if the update notification should be sent, or a false
   * if not.
   */
  template<class Trigger = TriggerAlways, class Counter = SimpleCounter>
  class ProgressCounter : public Counter {
  public:
    ProgressCounter(const Trigger& t) : trig(t) { }
    ProgressCounter(const Trigger& t, const Counter& c) : Counter(c), trig(t) { }

    //! Change the trigger
    void setTrigger(const Trigger& t) { trig = t; }

    //! Update the number of iterations and decide whether or not to
    //! notify based on the trigger policy
    void update(void) {
      ++Counter::count_;
      if (trig(this))
        Counter::notify();
    }

  private:
    Trigger trig;
  };



  // ------------------------------------------------------------------------------------------

  //! A basic progrss update "watcher", outputting dots for each
  //! update
  /**
   *
   * This class is the prototypic progress update class.  It has pre-
   * and post- condition output messages and then writes out a message
   * at each update.  You can control where the output goes by
   * instantiating with an appropriate ostream object.
   */
  class BasicProgress : public AbstractObserver {
  public:
    BasicProgress() : os_(std::cerr), prefix_("Progress - "), msg_("."), suffix_(" done!\n") { }
    BasicProgress(std::ostream& os, const std::string& prefix, const std::string& msg, const std::string& suffix) :
      os_(os), prefix_(prefix), msg_(msg), suffix_(suffix) { }

    BasicProgress(const std::string& prefix, const std::string& msg, const std::string& suffix) :
      os_(std::cerr), prefix_(prefix), msg_(msg), suffix_(suffix) { }

    virtual void start(SimpleCounter*) { os_ << prefix_; }
    virtual void update(SimpleCounter*) { os_ << msg_; }
    virtual void finish(SimpleCounter*) { os_ << suffix_; }

  protected:
    std::ostream& os_;
    std::string prefix_, msg_, suffix_;
  };


  //! Provide feedback by percent-complete with estimates of time
  //! remaining
  /**
   *
   * This class provides basically the same functionality as the
   * BasicProgress class, except that we override the update and
   * finish functions to provide more detailed information.  This
   * class requires a counter that implements the timeRemaining() and
   * fractionComplete() functions.
   */
  class PercentProgress : public BasicProgress {
  public:
    PercentProgress() : BasicProgress("Progress:\n", "complete", "") { }
    PercentProgress(std::ostream& os, const std::string& prefix, const std::string& msg, const std::string& suffix) :
      BasicProgress(os, prefix, msg, suffix) { }
    PercentProgress(const std::string& prefix, const std::string& msg, const std::string& suffix) :
      BasicProgress(std::cerr, prefix, msg, suffix) { }

    void update(SimpleCounter* s) {
      uint i = static_cast<uint>(floor(s->fractionComplete() * 100.0));
      os_ << i << "% " << msg_ << " (" << timeAsString(s->timeRemaining()) << " remaining)\n";
    }

    void finish(SimpleCounter* s) {
      os_ << suffix_;
      os_ << "Total elapsed time was " << timeAsString(s->elapsed()) << std::endl;
    }
  };



}



#endif
