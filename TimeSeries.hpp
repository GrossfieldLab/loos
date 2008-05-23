/*
  TimeSeries.hpp
  (c) 2008 Alan Grossfield


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Basic TimeSeries/vector math class
*/

#if !defined(TIMESERIES_HPP)
#define TIMESERIES_HPP

#include<vector>


using namespace std;

template<class T>
class tTimeSeries {
public:
    tTimeSeries() {
        init();
    }

    tTimeSeries(const vector<T> &inp) {
        _data = inp; 
    }

    tTimeSeries(const int size, const T*array) {
        _data.reserve(size);
        for (int i=0; i<size; i++) {
            _data.push_back(array[i]);
        }
    }

    tTimeSeries(const &tTimeSeries<T> inp) {
        _data = inp._data;
    }


    void init(void) {
        _data.clear();
    }

    void zero(void) {
        _data.assign(_data.size, (T)0.0);
    }

    T operator[](const int i) {
        return (_data.at(i)); // this way we get range checking
    }

    int size(void) const {
        return (_data.size());
    }

    T average(void) const {
        T ave = 0.0;
        for (int i=0; i < _data.size(); i++) {
            ave += _data[i];
        }
        ave /= _data.size();
        return (ave);
    }

    T variance(void) const {
        T ave = 0.0;
        T ave2 = 0.0;
        for (int i=0; i < _data.size(); i++) {
            ave += _data[i];
            ave2 += _data[i] * _data[i];
        }

        ave /= _data.size();
        ave2 /= _data.size();
        return (ave2 - ave*ave);
    }

    T stdev(void) const {
        return (sqrt(variance() ));
    }

    T sterr(void) const {
        // this assumes all points are statistically independent
        // Otherwise, needs to be multiplied by the square root of the 
        // correlation time, in units of the step interval for the time series
        return (stdev() / sqrt(_data.size()));  
    }

    tTimeSeries<T> running_average(void) const {
        // This is clumsy -- I'm copying the data 3 times, but it's
        // the only way I could think of to keep _data private, other than
        // making _data a pointer, which I may do at some point
        vector<T> tmp_data = _data;
        for (int i=1; i<_data.size(); i++) {
            _data[i] += _data[i-1];
            _data[i] /= (i+1);
        }

        tTimeSeries<T> new_series(*this);
        _data = tmp_data;
    }
    
private:
    vector<T> _data;
};

typedef tTimeSeries<double> TimeSeries;
typedef tTimeSeries<float> fTimeSeries;


#endif
