#include <iostream>
#include <string>
#include <stdexcept>

#include <loos_defs.hpp>
#include <timer.hpp>




#if !defined(PROGRESSCOUNTERS_HPP)
#define PROGRESSCOUNTERS_HPP

namespace loos {

  class SimpleCounter;


  class AbstractObserver {
  public:
    virtual ~AbstractObserver() { }
    
    virtual void start(SimpleCounter*) =0;
    virtual void finish(SimpleCounter*) =0;
    virtual void update(SimpleCounter*) =0;
  };


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

    virtual void notify(void) {
      Observers::iterator i;
      for (i = observers.begin(); i != observers.end(); ++i)
        (*i)->update(this);
    }

    virtual void finish(void) {
      timer_.stop();
      Observers::iterator i;
      for (i = observers.begin(); i != observers.end(); ++i)
        (*i)->finish(this);
    }


    virtual void start(void) {
      count_ = 0;
      timer_.start();
      Observers::iterator i;
      for (i = observers.begin(); i != observers.end(); ++i)
        (*i)->start(this);
    }


    uint count(void) const { return(count_); }

    virtual double elapsed(void) { return(timer_.elapsed()); }

    virtual uint remaining(void) { throw(std::logic_error("remaining() is unimplemented")); }
    virtual double timeRemaning(void) { throw(std::logic_error("timeRemaining() is unimplemented")); }
    virtual double fractionComplete(void) { throw(std::logic_error("fractionComplete() is unimplemented")); }

  protected:
    uint count_;
    Timer<WallTimer> timer_;
    Observers observers;
  };



  class EstimatingCounter : public SimpleCounter {
  public:
    EstimatingCounter(const uint n) : expected(n) { }
    void setExpected(uint n) { expected = n; }

    uint remaining(void) { return(expected - count_); }

    double timeRemaining(void) {
      double avg = timer_.elapsed() / count_;
      return(remaining() * avg);
    }

    double fractionComplete(void) { return(static_cast<double>(count_) / expected); }

  protected:
    uint expected;
  };


  class TriggerAlways {
  public:
    bool operator()(SimpleCounter* s) { return(true); }
  };



  template<class Trigger = TriggerAlways, class Counter = SimpleCounter>
  class ProgressCounter : public Counter {
  public:
    ProgressCounter(Trigger& t) : trig(t) { }
    ProgressCounter(Trigger& t, const Counter& c) : Counter(c), trig(t) { }

    void setTrigger(const Trigger& t) { trig = t; }

    void update(void) {
      ++Counter::count_;
      if (trig(this))
        Counter::notify();
    }

  private:
    Trigger trig;
  };



  // ------------------------------------------------------------------------------------------

  class DotProgress : public AbstractObserver {
  public:
    DotProgress(std::ostream& os, const std::string& prefix, const std::string& msg, const std::string& suffix) :
      os_(os), prefix_(prefix), msg_(msg), suffix_(suffix) { }

    DotProgress(const std::string& prefix, const std::string& msg, const std::string& suffix) :
      os_(std::cerr), prefix_(prefix), msg_(msg), suffix_(suffix) { }

    void start(SimpleCounter* subj) { os_ << prefix_; }
    void update(SimpleCounter* subj) { os_ << msg_; }
    void finish(SimpleCounter* subj) { os_ << suffix_; }

  private:
    std::ostream& os_;
    std::string prefix_, msg_, suffix_;
  };


  class PercentProgress : public AbstractObserver {
  public:
    PercentProgress(std::ostream& os, const std::string& prefix, const std::string& msg, const std::string& suffix) :
      os_(os), prefix_(prefix), msg_(msg), suffix_(suffix) { }

    PercentProgress(const std::string& prefix, const std::string& msg, const std::string& suffix) :
      os_(std::cerr), prefix_(prefix), msg_(msg), suffix_(suffix) { }

    void start(SimpleCounter* s) { os_ << prefix_; }
    void update(SimpleCounter* s) {
      uint i = static_cast<uint>(floor(s->fractionComplete() * 100.0));
      os_ << i << "% " << msg_ << "(" << timeAsString(s->elapsed()) << " remaining)\n";
    }

    void finish(SimpleCounter* s) {
      os_ << suffix_;
      os_ << "Total elapsed time was " << timeAsString(s->elapsed()) << std::endl;
    }
  private:
    std::ostream& os_;
    std::string prefix_, msg_, suffix_;
  };



}



#endif
