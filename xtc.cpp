#include <xtc.hpp>



namespace loos {

  namespace internal {

    namespace xtc {


      int sizeofint(int size) {
        int n = 0;
        while ( size > 0 ) {
          size >>= 1;
          ++n;
        }
        return(n);
      }


      int sizeofints(uint* data, const uint n) {
        uint nbytes = 1;
        uint bytes[32];
        uint nbits = 0;
    
        bytes[0] = 1;
        for (uint i = 0; i < n; ++i) {
          uint tmp = 0;
          uint bytecnt;
          for (bytecnt = 0; bytecnt < nbytes; ++bytecnt) {
            tmp += bytes[bytecnt] * data[i];
            bytes[bytecnt] = tmp & 0xff;
            tmp >>= 8;
          }
          while (tmp != 0) {
            bytes[bytecnt++] = tmp & 0xff;
            tmp >>= 8;
          }
          nbytes = bytecnt;
        }

        uint num = 1;
        --nbytes;
        while (bytes[nbytes] >= num) {
          ++nbits;
          num *= 2;
        }

        return(nbits + nbytes*8);
      }


      int decodebits(int* buf, uint nbits) {

        int mask = (1 << nbits) -1;

        unsigned char *cbuf = reinterpret_cast<unsigned char*>(buf) + 3*sizeof(*buf);
        int cnt = buf[0];
        uint lastbits = static_cast<uint>(buf[1]);
        uint lastbyte = static_cast<uint>(buf[2]);
    
        int num = 0;
        while (nbits >= 8) {
          lastbyte = ( lastbyte << 8 ) | cbuf[cnt++];
          num |=  (lastbyte >> lastbits) << (nbits - 8);
          nbits -=8;
        }
        if (nbits > 0) {
          if (lastbits < nbits) {
            lastbits += 8;
            lastbyte = (lastbyte << 8) | cbuf[cnt++];
          }
          lastbits -= nbits;
          num |= (lastbyte >> lastbits) & ((1 << nbits) -1);
        }
        num &= mask;
        buf[0] = cnt;
        buf[1] = static_cast<int>(lastbits);
        buf[2] = static_cast<int>(lastbyte);

        return(num);
      }

      void decodeints(int* buf, const int num_of_ints, int num_of_bits,
                      uint* sizes, int* nums) {
        int bytes[32];
        int i, j, num_of_bytes, p, num;
  
        bytes[1] = bytes[2] = bytes[3] = 0;
        num_of_bytes = 0;
        while (num_of_bits > 8) {
          bytes[num_of_bytes++] = decodebits(buf, 8);
          num_of_bits -= 8;
        }
        if (num_of_bits > 0) {
          bytes[num_of_bytes++] = decodebits(buf, num_of_bits);
        }
        for (i = num_of_ints-1; i > 0; i--) {
          num = 0;
          for (j = num_of_bytes-1; j >=0; j--) {
            num = (num << 8) | bytes[j];
            p = num / sizes[i];
            bytes[j] = p;
            num = num - p * sizes[i];
          }
          nums[i] = num;
        }
        nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
      }
    }

  }


}
