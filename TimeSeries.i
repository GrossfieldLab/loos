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



    void init(void) {
      _data.clear();
    }

    void zero(void) {
      _data.assign(_data.size(), (T)0.0);
    }


    unsigned int size(void) const {
      return (_data.size());
    }


    TimeSeries<T> operator+=(const T val) {
      for (unsigned int i=0; i<_data.size(); i++) {
        _data[i] += val;
      }
      return(*this);
    }

    TimeSeries<T> operator+=(const TimeSeries<T> &rhs) {
      if (_data.size() != rhs.size())
        throw(std::runtime_error("mismatched timeseries sizes in +="));
      for (unsigned int i=0; i<_data.size(); i++) {
        _data[i] += rhs[i];
      }
      return(*this);
    }

    TimeSeries<T> operator+(const T val) const {
      TimeSeries<T> res(*this);
      for (unsigned int i=0; i<_data.size(); i++) {
        res[i] += val;
      }
      return(res);
    }


    // friend TimeSeries<T> operator+(const T lhs, const TimeSeries<T> &rhs) {
    //   TimeSeries<T> res( rhs.size(), (T) 0.0 );
    //   for (unsigned int i=0; i<rhs.size(); i++) {
    //     res[i] = lhs + rhs[i];
    //   }
    //   return(res);
    // }

    TimeSeries<T> operator+(const TimeSeries<T> &rhs) const {
      if (_data.size() != rhs.size())
        throw(std::runtime_error("mismatched timeseries sizes in +="));

      TimeSeries<T> res(*this);
      for (unsigned int i=0; i<res.size(); i++) {
        res[i] += rhs[i];
      }
      return(res);
    }

    TimeSeries<T> operator-=(const T val) {
      for (unsigned int i=0; i<_data.size(); i++) {
        _data[i] -= val;
      }
      return(*this);
    }

    TimeSeries<T> operator-=(const TimeSeries<T> &rhs) {
      if (_data.size() != rhs.size())
        throw(std::runtime_error("mismatched sizes of time series"));
      for (unsigned int i=0; i<_data.size(); i++) {
        _data[i] -= rhs[i];
      }
      return(*this);
    }

    TimeSeries<T> operator-(const T val) const {
      TimeSeries<T> res(*this);
      for (unsigned int i=0; i<_data.size(); i++) {
        res[i] -= val;
      }
      return(res);
    }

    TimeSeries<T> operator-(const TimeSeries<T> &rhs) const {
      if (_data.size() != rhs.size())
        throw(std::runtime_error("mismatched timeseries sizes in +="));

      TimeSeries<T> res(*this);
      for (unsigned int i=0; i<res.size(); i++) {
        res[i] -= rhs[i];
      }
      return(res);
    }



    TimeSeries<T> operator-() const {
      TimeSeries<T> res(*this);
      for (unsigned int i=0; i<_data.size(); i++) {
        res[i] = -res[i];
      }
      return(res);
    }

    // friend TimeSeries<T> operator-(const T lhs, const TimeSeries<T> &rhs) {
    //   TimeSeries<T> res( rhs.size() );
    //   for (unsigned int i=0; i<rhs.size(); i++) {
    //     res[i] = lhs - rhs[i];
    //   }
    //   return(res);
    // }

    TimeSeries<T> operator*=(const T val) {
      for (unsigned int i=0; i<_data.size(); i++) {
        _data[i] *= val;
      }
      return(*this);
    }

    TimeSeries<T> operator*(const T val) const {
      TimeSeries<T> res(*this);
      for (unsigned int i=0; i<_data.size(); i++) {
        res[i] *= val;
      }
      return(res);
    }

    TimeSeries<T> operator*=(const TimeSeries<T> &rhs) {
      if ( _data.size() != rhs.size() ) 
        throw(std::runtime_error("mismatched timeseries sizes in *="));
      for (unsigned int i=0; i<_data.size(); i++) {
        _data[i] *= rhs[i];
      }
      return(*this);
    }

    TimeSeries<T> operator*(const TimeSeries<T> &rhs) const {
      if ( _data.size() != rhs.size() ) 
        throw(std::runtime_error("mismatched timeseries sizes in *="));

      TimeSeries<T> res(*this);
      for (unsigned int i=0; i<_data.size(); i++) {
        res[i] *= rhs[i];
      }
      return(res);
    }

    // friend TimeSeries<T> operator*(const T lhs, const TimeSeries<T> &rhs) {
    //   TimeSeries<T> res( rhs.size(), (T) 0.0 );
    //   for (unsigned int i=0; i<rhs.size(); i++) {
    //     res[i] = lhs * rhs[i];
    //   }
    //   return(res);
    // }


    TimeSeries<T> operator/=(const T val) {
      for (unsigned int i=0; i<_data.size(); i++) {
        _data[i] /= val;
      }
      return(*this);
    }

    TimeSeries<T> operator/(const T val) const {
      TimeSeries<T> res(*this);
      for (unsigned int i=0; i<_data.size(); i++) {
        res[i] /= val;
      }
      return(res);
    }

    TimeSeries<T> operator/=(const TimeSeries<T> &rhs) {
      if ( _data.size() != rhs.size() ) 
        throw(std::runtime_error("mismatched timeseries sizes in *="));
      for (unsigned int i=0; i<_data.size(); i++) {
        _data[i] /= rhs[i];
      }
      return(*this);
    }

    TimeSeries<T> operator/(const TimeSeries<T> &rhs) const {
      if ( _data.size() != rhs.size() ) 
        throw(std::runtime_error("mismatched timeseries sizes in *="));

      TimeSeries<T> res(*this);
      for (unsigned int i=0; i<_data.size(); i++) {
        res[i] /= rhs[i];
      }
      return(res);
    }

    // friend TimeSeries<T> operator/(const T lhs, const TimeSeries<T> &rhs) {
    //   TimeSeries<T> res( rhs.size() );
    //   for (unsigned int i=0; i<rhs.size(); i++) {
    //     res[i] = lhs / rhs[i];
    //   }
    //   return(res);
    // }

    TimeSeries<T> copy(void) const {
      TimeSeries<T> result(_data.size(), 0.0);
      for (unsigned int i = 0; i < _data.size(); i++) {
        result[i] = _data[i];
      }
      return(result);
    }

    //! Remove num_points from the front of the time series, as you
    //! would to eliminate the equilibration time
    void set_skip(unsigned int num_points) {
        if (num_points < _data.size()) {
            _data.erase(_data.begin(), _data.begin()+num_points);
        } 
        else {
            throw(std::out_of_range("num_points has invalid value in set_skip"));
        }
    }

    //! Return average of time series
    T average(void) const {
      T ave = 0.0;
      for (unsigned int i=0; i < _data.size(); i++) {
        ave += _data[i];
      }
      ave /= _data.size();
      return (ave);
    }

    //! Return variance of time series.
    T variance(void) const {
      T ave = 0.0;
      T ave2 = 0.0;
      for (unsigned int i=0; i < _data.size(); i++) {
        ave += _data[i];
        ave2 += _data[i] * _data[i];
      }

      ave /= _data.size();
      ave2 /= _data.size();
      return (ave2 - ave*ave);
    }

    //! Return standard deviation of time series.
    T stdev(void) const {
      return (sqrt(variance() ));
    }

    //! Return standard error of time series.
    //! This assumes all points are statistically independent
    //! Otherwise, needs to be multiplied by the square root of the 
    //! correlation time, in units of the step interval for the time series
    T sterr(void) const {
      return (stdev() / sqrt(_data.size()));  
    }

    //! Return a new timeseries of the same size as the current one,
    //! containing the running average of the time series
    TimeSeries<T> running_average(void) const {
      TimeSeries<T> result(_data.size());
      T sum = 0.0;
      for (unsigned int i=0; i<_data.size(); i++) {
        sum += _data[i];
        result[i] = sum / (i+1);
      }
      return(result);

    }

    //! Return a new timeseries containing the windowed average.
    //! ith value of the new time series =  1/window * sum(data[i:i+window]).
    //! NOTE: The present algorithm is relatively fast, but can be prone
    //! to roundoff.
    TimeSeries<T> windowed_average(const uint window) const {

      if (window > _data.size() )
        throw(std::out_of_range("Error in windowed_average: window too large"));

      TimeSeries<T> result(_data.size() - window);
      T sum = 0;
      for (uint i=0; i<window; i++) {
        sum += _data[i];
      }
      result[0] = sum / window;


      for (unsigned int i=1; i < result.size(); i++) {
        sum = sum - _data[i-1] + _data[i+window-1];
        result[i] = sum / window;
      }

      return(result);
    }

    //! Return the variance of the block average for the time series.
    //! Divides the timeseries into num_blocks equally sized blocks
    //! (discarding the remaining blocks at the end), computes the average
    //! for each block, and returns the variance of the averages.
    //! This is useful for doing Flyvjberg and Petersen-style block averaging.
    //! Flyvbjerg, H. & Petersen, H. G. J. Chem. Phys., 1989, 91, 461-466
    // 
    T block_var(const int num_blocks) const {
      int points_per_block = size() / num_blocks;
      T block_ave = 0.0;
      T block_ave2 = 0.0;
      for (int i=0; i<num_blocks; i++) {
        T block_sum = 0.0;
        int offset = i*points_per_block;
        for (int j=0; j< points_per_block; j++) {
          block_sum += _data[offset + j];
        }
        T ave = block_sum / points_per_block;
        block_ave += ave;
        block_ave2 += ave*ave;
      }
      block_ave /= num_blocks;
      block_ave2 /= num_blocks;

      // The variance must be computed with N-1, not N, because
      // we determine the mean from the data (as opposed to independently
      // specifying it)
      T ratio = num_blocks/(num_blocks-1.0);
      return (block_ave2 - block_ave*block_ave)*ratio;
    }

    TimeSeries<T> correl(const int max_time, 
                         const int interval=1, 
                         const bool normalize=true,
                         T tol=1.0e-8) const {

      TimeSeries<T> data = copy();
      uint n = abs(max_time);
      if (n > data.size()) {
        throw(std::runtime_error("Can't take correlation time longer than time series"));
      }

      n /= interval;
      TimeSeries<T> c(n, 0.0);

      // normalize the data
      if (normalize) {
          data -= data.average();
          T dev = data.stdev();
            
          // drop through if this is a constant array
          if (dev < tol) {
            c._data.assign(n, 1.0);
            return(c);
          }
            
          data /= dev;
      }

      std::vector<int> num_pairs(n);
      num_pairs.assign(n, 0);
      // TODO: This is the O(N^2) way -- there are much faster
      // algorithms for long time series
      for (int i = 0; i < max_time; i+=interval) {
        int index = i / interval;
        for (unsigned int j = 0; j < data.size() - i; j++) {
          c[index] += data[j] * data[j+i];
          num_pairs[index]++;
        }

      }

      // Divide each value by the number of pairs used to generated it
      for (uint i = 0; i < n; i++) {
        c[i] /= num_pairs[i];
      }

      return(c);
    }

  // Vector interface...
  void push_back(const T& x) { _data.push_back(x); }



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


  %template(DTimeSeries) loos::TimeSeries<double>;
  



}
