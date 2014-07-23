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

    XTCWriter() :
      buf1size(0), buf2size(0),
      buf1(0), buf2(0)
    {
      
    }
    

    ~XTCWriter() {
      if (buf1)
	delete[] buf1;
      if (buf2)
	delete[] buf2;
    }


  private:
    int sizeofint(const int size) const;
    int sizeofints(const int num_of_bits, const unsigned int sizes[]) const;
    void encodebits(int* buf, int num_of_bits, const int num) const;
    void encodeints(int* buf, const int num_of_ints, const int num_of_bits,
		    const unsigned int* sizes, const unsigned int* nums) const;
    int writeCompressedCoordsFloat(float* ptr, int size, float precision);
    int writeCompressedCoordsDouble(double* ptr, int size, double precision);
       
    void allocateBuffers(const size_t size);
    

  private:
    uint buf1size, buf2size;
    float* buf1;
    float* buf2;
  };


};



#endif
