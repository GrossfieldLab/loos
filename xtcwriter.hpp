// (c) 2014 Tod D. Romo, Grossfield Lab, URMC

#if !defined(LOOS_XTCWRITER_HPP)
#define LOOS_XTC_WRITER_HPP

#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>
#include <xdr.hpp>

namespace loos {

  class XTCWriter {
    static const int magicints[];
    static const int firstidx;
    
    static const int lastidx;
    

  public:


  private:
    int sizeofint(const int size) const;
    int sizeofints(const int num_of_bits, const unsigned int sizes[]) const;
    void encodebits(int* buf, int num_of_bits, const int num) const;
    void encodeints(int* buf, const int num_of_ints, const int num_of_bits,
		    const unsigned int* sizes, const unsigned int* nums) const;
    int writeCompressedCoordsFloat(float* ptr, int size, float precision);
    int writeCompressedCoordsDouble(double* ptr, int size, double precision);
    

  private:
    uint buf1size, buf2size;
    float* buf1;
    float* buf2;
  };


};



#endif
