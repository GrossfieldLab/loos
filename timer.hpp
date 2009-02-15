// A simple timer class...

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008,2009, Tod D. Romo, Alan Grossfield
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


#if !defined(TIMER_HPP)
#define TIMER_HPP

#include <iostream>
#include <string>
#include <stdexcept>

#include <sys/time.h>
#include <sys/resource.h>
#include <boost/format.hpp>

#include <loos_defs.hpp>


namespace loos {

  //! Policy class for tracking wall-time
  class WallTimer {
  public:
    double currentTime(void) const {
      struct timeval tv;
      
      int i = gettimeofday(&tv, 0);
      if (i < 0)
        throw(std::runtime_error("Error in gettimeofday()"));
      return(static_cast<double>(tv.tv_sec) + static_cast<double>(tv.tv_usec) * 1e-6);
    }
  };
  
  //! Policy class for tracking only user process time
  class UserTimer {
  public:
    double currentTime(void) const {
      struct rusage ru;
      int i = getrusage(RUSAGE_SELF, &ru);
      if (i < 0)
        throw(std::runtime_error("Error in getrusage()"));
      return(static_cast<double>(ru.ru_utime.tv_sec) + static_cast<double>(ru.ru_utime.tv_usec) * 1e-6);
    }
  };


  //! Class for tracking time
  /**
   * Can either track wall-time or user-time, depending on the timer
   * policy class used...  The default is to track wall-time.
   * Starting and stopping the timer does not actually do anything
   * other than adjust how the current time is tracked, i.e. no real
   * timers (interrupts, etc) are created.
   *
   * Time is tracked to microsecond precision...  See gettimeofday(2)
   * or getrusage(2) for more information...
   */

  template<class TimerType = WallTimer>
  class Timer : public TimerType {
  public:
    Timer() : t0(0), t1(0), avg(0), lapt(0), n(0), running(false) { }

    //! Starts the timer
    void start(void) {
      lapt = t0 = TimerType::currentTime();
      n=0;
      avg = 0.0;
      running = true;
    }

    //! Stops the timer
    // Automatically adds a lap...
    double stop(void) {
      t1 = TimerType::currentTime();
      avg += (t1 - lapt);
      ++n;
      running = false;
      return(t1-t0);
    }

    //! Return the elapsed time.
    /**
     * If the timer is running, it returns the difference between
     * current time and when the timer was started.  If the timer is
     * stopped, then it returns the length of time the timer was
     * active.
     */
    double elapsed(void) const { return( running ? TimerType::currentTime() - t0 : t1 - t0 ); }

    //! Returns the current lap time
    /**
     * Lap time is a mechanism for tracking interval times.  Each time
     * lap() is called, it tracks not only the current lap-time but
     * the average time for each lap.  When the timer is stopped, a
     * lap is automatically added.
     */
    double lap(void) {
      if (!running)    // Don't adjust lap time if the timer ain't running...
        return(0.0);

      double t = TimerType::currentTime();
      double lt = t - lapt;
      lapt = t;
      avg += lt;
      ++n;
    }

    //! Return the current average lap-time...
    double averageLapTime(void) { return(avg / n); }
    
  private:
    double t0, t1, avg, lapt;
    ulong n;
    bool running;
  };

  template<class T>
  std::ostream& operator<<(std::ostream& os, const Timer<T>& t) {
    std::string s = timeAsString(t.elapsed());
    os << "Elapsed time " << s;
    return(os);
  }
}




#endif
