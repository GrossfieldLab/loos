#if !defined(XTC_HPP)
#define XTC_HPP

#include <loos_defs.hpp>

#include <xdr.hpp>
#include <Trajectory.hpp>

#include <boost/format.hpp>

namespace loos {

  namespace internal {
    namespace xdr {

      int sizeofint(int);
      int sizeofints(uint*, const int);
      int decodebits(int*, int);
      void decodeints(int*, const int, int, uint*, int*);
    
  
      const int magicints[] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
        80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024, 1290, // 31
        1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003, // 41
        16384, 20642, 26007, 32768, 41285, 52015, 65536,82570, 104031, // 50
        131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561, // 58
        832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021, // 65
        4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216 // 72
      };
      
      const int firstidx = 9;
      const int lastidx = 73;
      
    }

  }
  
  using namespace internal::xdr;


  template<typename T>
  int xdrfile_read_compr_coord(std::vector<T> fp, /* length 3*ncoord */ 
                                     float     *precision,
                                     internal::XDR*   xfp)
   {
     int minint[3], maxint[3], *lip;
     int smallidx, minidx, maxidx;
     uint sizeint[3], sizesmall[3], bitsizeint[3], size3;
     int k, *buf1, *buf2, lsize, flag;
     int smallnum, smaller, larger, i, is_smaller, run;
     T inv_precision;
     int tmp, *thiscoord,  prevcoord[3];
     unsigned int bitsize;
  
     
     if (!xfp->read(lsize))
       return(0);

     size3 = lsize * 3;

     /* Dont bother with compression for three atoms or less */
     if(lsize<=9) {
       T* tmp = new T[size3];
       xfp->read(tmp, size3);
       for (int i=0; i<size3; ++i)
         fp.push_back(tmp[i]);
       delete[] tmp;
       return(lsize);
     }

     /* Compression-time if we got here. Read precision first */
     xfp->read(precision);
  
     int size3padded = size3 * 1.2;
     buf1 = new int[size3padded];
     buf2 = new int[size3padded];
     /* buf2[0-2] are special and do not contain actual data */
     buf2[0] = buf2[1] = buf2[2] = 0;
     xfp->read(minint, 3);
     xfp->read(maxint, 3);
  
     sizeint[0] = maxint[0] - minint[0]+1;
     sizeint[1] = maxint[1] - minint[1]+1;
     sizeint[2] = maxint[2] - minint[2]+1;
	
     /* check if one of the sizes is to big to be multiplied */
     if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff) {
       bitsizeint[0] = sizeofint(sizeint[0]);
       bitsizeint[1] = sizeofint(sizeint[1]);
       bitsizeint[2] = sizeofint(sizeint[2]);
       bitsize = 0; /* flag the use of large sizes */
     } else {
       bitsize = sizeofints(sizeint, 3);
     }
	
     if (!xfp->read(smallidx)) {
       delete[] buf1;
       delete[] buf2;
       return(0);
     }

     tmp=smallidx+8;
     maxidx = (lastidx<tmp) ? lastidx : tmp;
     minidx = maxidx - 8; /* often this equal smallidx */
     tmp = smallidx-1;
     tmp = (firstidx>tmp) ? firstidx : tmp;
     smaller = magicints[tmp] / 2;
     smallnum = magicints[smallidx] / 2;
     sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
     larger = magicints[maxidx];

     /* buf2[0] holds the length in bytes */
  
     if (!xfp->read(buf2, 1)) {
       delete[] buf1;
       delete[] buf2;
       return(0);
     }

     if (!xfp->read(reinterpret_cast<char*>(&(buf2[3])), static_cast<uint>(buf2[0]))) {
       delete[] buf1;
       delete[] buf2;
       return(0);
     }
     buf2[0] = buf2[1] = buf2[2] = 0;
  
     inv_precision = 1.0 / * precision;
     run = 0;
     i = 0;
     lip = buf1;
     while ( i < lsize ) {
       thiscoord = (int *)(lip) + i * 3;
    
       if (bitsize == 0) {
         thiscoord[0] = decodebits(buf2, bitsizeint[0]);
         thiscoord[1] = decodebits(buf2, bitsizeint[1]);
         thiscoord[2] = decodebits(buf2, bitsizeint[2]);
       } else {
         decodeints(buf2, 3, bitsize, sizeint, thiscoord);
       }
    
       i++;
       thiscoord[0] += minint[0];
       thiscoord[1] += minint[1];
       thiscoord[2] += minint[2];
    
       prevcoord[0] = thiscoord[0];
       prevcoord[1] = thiscoord[1];
       prevcoord[2] = thiscoord[2];
    
       flag = decodebits(buf2, 1);
       is_smaller = 0;
       if (flag == 1) {
         run = decodebits(buf2, 5);
         is_smaller = run % 3;
         run -= is_smaller;
         is_smaller--;
       }
       if (run > 0) {
         thiscoord += 3;
         for (k = 0; k < run; k+=3) {
           decodeints(buf2, 3, smallidx, sizesmall, thiscoord);
           i++;
           thiscoord[0] += prevcoord[0] - smallnum;
           thiscoord[1] += prevcoord[1] - smallnum;
           thiscoord[2] += prevcoord[2] - smallnum;
           if (k == 0) {
             /* interchange first with second atom for better
              * compression of water molecules
              */
             tmp = thiscoord[0]; thiscoord[0] = prevcoord[0];
             prevcoord[0] = tmp;
             tmp = thiscoord[1]; thiscoord[1] = prevcoord[1];
             prevcoord[1] = tmp;
             tmp = thiscoord[2]; thiscoord[2] = prevcoord[2];
             prevcoord[2] = tmp;

             T x = prevcoord[0] * inv_precision;
             T y = prevcoord[1] * inv_precision;
             T z = prevcoord[2] * inv_precision;
             fp.push_back(x);
             fp.push_back(y);
             fp.push_back(z);
             std::cerr << "(" << x << "," << y << "," << z << ")\n";
           } else {
             prevcoord[0] = thiscoord[0];
             prevcoord[1] = thiscoord[1];
             prevcoord[2] = thiscoord[2];
           }
           T x = thiscoord[0] * inv_precision;
           T y = thiscoord[1] * inv_precision;
           T z = thiscoord[2] * inv_precision;
           fp.push_back(x);
           fp.push_back(y);
           fp.push_back(z);
           std::cerr << "(" << x << "," << y << "," << z << ")\n";

         }
       } else {
         T x = thiscoord[0] * inv_precision;
         T y = thiscoord[1] * inv_precision;
         T z = thiscoord[2] * inv_precision;
         fp.push_back(x);
         fp.push_back(y);
         fp.push_back(z);
         std::cerr << "(" << x << "," << y << "," << z << ")\n";
       }
       smallidx += is_smaller;
       if (is_smaller < 0) {
         smallnum = smaller;
         if (smallidx > firstidx) {
           smaller = magicints[smallidx - 1] /2;
         } else {
           smaller = 0;
         }
       } else if (is_smaller > 0) {
         smaller = smallnum;
         smallnum = magicints[smallidx] / 2;
       }
       sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
     }

     delete[] buf1;
     delete[] buf2;
     return(lsize);
   }



};


#endif
