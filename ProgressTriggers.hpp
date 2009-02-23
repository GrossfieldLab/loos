#include <loos_defs.hpp>
#include <ProgressCounters.hpp>

#if !defined(PROGRESSTRIGGERS_HPP)
#define PROGRESSTRIGGERS_HPP

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
