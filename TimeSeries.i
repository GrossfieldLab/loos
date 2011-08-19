/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Alan Grossfield
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


%header %{
#include <TimeSeries.hpp>
%}

namespace loos {

  //! Time Series Class
  /*!
   *  This class provides basic operations on a time series, such
   *  as averaging, standard deviation, etc
   *
   *  One can do standard arimethic operations on time series combined with
   *  scalars or other timeseries (as long as the two time series are the same
   *  length).
   *
   */

template<class T>
class TimeSeries {
public:
  typedef typename std::vector<T>::iterator       iterator;
  typedef typename std::vector<T>::const_iterator const_iterator;
  typedef const T&  const_reference;
  typedef T&        reference;


    TimeSeries() {
      init();
    }

    TimeSeries(const std::vector<T> &inp) {
      _data = inp; 
    }

    TimeSeries(const uint size, const T *array) {
      _data.reserve(size);
      for (unsigned int i=0; i<size; i++) {
        _data.push_back(array[i]);
      }
    }

    TimeSeries(const TimeSeries<T> &inp) {
      _data = inp._data;
    }

    TimeSeries(const uint n) {
      _data = std::vector<T>(n, 0);
    }

    TimeSeries(const uint n, const T val) {
      _data.assign(n, (T) val);
    }

    //! Resize the TimeSeries by calling the underlying vector's resize
    void resize(const uint n, const T val= (T) 0.0) {
        _data.resize(n, val);
    }


    //! Read a simple text file and create a timeseries
    //! The file is assumed to be simple columnated data.  Blank lines and 
    //! lines starting with "#" are ignored.
    TimeSeries (const std::string &filename, const int col=2) {
        std::ifstream ifs(filename.c_str());
        if (!ifs) {
            throw(std::runtime_error("Cannot open timeseries file " 
                                     + filename));
        }

        std::string line;
        while (ifs.good()) {
            getline(ifs, line);
            if ( (line.substr(0,1) == "#") || (line.empty()) ){
                // comment -- do nothing
            }
            else {
                std::stringstream s(line);
                double val;
                int i=0;
                while (i < col) {
                    if (!s.good()) {
                  throw(std::runtime_error("Problem reading timeseries file "
                                           + filename));
                    }
                    s >> val;
                    i++;
                }
                _data.push_back(val);    
            }
        }
    }



  void zero(void);
  unsigned int size(void) const;
  TimeSeries<T> operator+=(const T val);

  TimeSeries<T> operator+=(const TimeSeries<T> &rhs);
  TimeSeries<T> operator+(const T val) const;
  TimeSeries<T> operator+(const TimeSeries<T> &rhs) const;
  TimeSeries<T> operator-=(const T val);
  TimeSeries<T> operator-=(const TimeSeries<T> &rhs);
  TimeSeries<T> operator-(const T val) const;
  TimeSeries<T> operator-(const TimeSeries<T> &rhs) const;
  TimeSeries<T> operator-() const;
  TimeSeries<T> operator*=(const T val);
  TimeSeries<T> operator*(const T val) const;
  TimeSeries<T> operator*=(const TimeSeries<T> &rhs);
  TimeSeries<T> operator*(const TimeSeries<T> &rhs) const;
  TimeSeries<T> operator/=(const T val);
  TimeSeries<T> operator/(const T val) const;
  TimeSeries<T> operator/=(const TimeSeries<T> &rhs);
  TimeSeries<T> operator/(const TimeSeries<T> &rhs) const;
  TimeSeries<T> copy(void) const;
  void set_skip(unsigned int num_points);
  T average(void) const;
  T variance(void) const;
  T stdev(void) const;
  T sterr(void) const;
  TimeSeries<T> running_average(void) const;
  TimeSeries<T> windowed_average(const uint window) const;
  T block_var(const int num_blocks) const;
  TimeSeries<T> correl(const int max_time, 
                       const int interval=1, 
                       const bool normalize=true,
                       T tol=1.0e-8) const;
  void push_back(const T& x);



};


  %extend TimeSeries<double> {

    double __getitem__(const int i) {
      return((*$self)[i]);
    }
    
    void __setitem__(const int i, const double d) {
      (*$self)[i] = d;
    }
  };

  %rename(__add__)  loos::TimeSeries<double>::operator+;
  %rename(__sub__) loos::TimeSeries<double>::operator-;
  %rename(__mul__) loos::TimeSeries<double>::operator*;
  %rename(__div__) loos::TimeSeries<double>::operator/;


  %template(TimeSeriesDbl) loos::TimeSeries<double>;
  



}
