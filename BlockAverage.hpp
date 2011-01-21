#if !defined(BLOCKAVERAGE_HPP)
#define BLOCKAVERAGE_HPP

#include "loos_defs.hpp"
#include <vector>

#include "exceptions.hpp"


namespace loos {

  template<typename F>
  double blockStandardError(const F& op, const uint start, const uint end, const uint blocksize) {
    uint nblocks = (end - start + 1) / blocksize;

    std::vector<double> values;
    double mean = 0.0;
    for (uint block = 0; block < nblocks; ++block) {
      double v = op(start + block * blocksize, blocksize);
      mean += v;
      values.push_back(v);
    }
    mean /= nblocks;

    double var = 0.0;
    for (std::vector<double>::iterator i = values.begin(); i != values.end(); ++i) {
      double d = *i - mean;
      var += d*d;
    }
    var /= (nblocks - 1);

    return(sqrt(var / nblocks));
  }


  template<typename F>
  std::vector<double> blockAverage(const F& op, const uint start, const uint end, const std::vector<uint>& block_sizes) {

    std::vector<double> errors;
    for (std::vector<uint>::const_iterator i = block_sizes.begin(); i != block_sizes.end(); ++i)
      errors.push_back(blockStandardError(op, start, end, *i));

    return(errors);
  }


  template<typename T>
  struct VectorBlockAverage {
    VectorBlockAverage(const std::vector<T>& d) : data(d) { }

    double operator()(const uint start, const uint blocksize) {
      if (start < 0 || start+blocksize >= data.size())
        throw(LOOSError("Invalid parameters to VectorBlockAverage::operator()"));

      double mean = 0.0;
      for (uint i=start; i < start+blocksize; ++i)
        mean += data[i];
      return(mean / blocksize);
    }


    const std::vector<T>& data;
  };

  template<typename T>
  std::vector<double> blockAverage(const std::vector<T>& data, const std::vector<uint>& block_sizes) {
    VectorBlockAverage<T> V(data);

    return(blockAverage(V, 0, data.size(), block_sizes));
  }


}



#endif
