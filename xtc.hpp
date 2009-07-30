#if !defined(XTC_HPP)
#define XTC_HPP
#include <ios>
#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

#include <loos_defs.hpp>

#include <xdr.hpp>
#include <Trajectory.hpp>

#include <boost/format.hpp>

namespace loos {

  template<typename T>
  class XTC : public Trajectory {

    struct Header {
      int magic;
      int natoms, step;
      float time, box[9];
    };

  public:
    explicit XTC(const std::string& s) : Trajectory(s), xdr_file(ifs()),natoms_(0) { scanFrames(); }
    explicit XTC(const char* p) : Trajectory(p), xdr_file(ifs()),natoms_(0) { scanFrames(); } 

    uint natoms(void) const { return(natoms_); }
    float timestep(void) const { return(0); }
    uint nframes(void) const { return(frame_indices.size()); }
    bool hasPeriodicBox(void) const { return(true); }
    GCoord periodicBox(void) const { return(box); }
    std::vector<GCoord> coords(void) { std::vector<GCoord> v; return(v); }
    void updateGroupCoords(AtomicGroup& g) { }
    bool parseFrame(void) { return(false); }

    void seekNextFrameImpl(void) { }
    void seekFrameImpl(uint i);
    void rewindImpl(void) { }

  private:
    static const int magicints[];
    static const int firstidx;
    static const int lastidx;

    internal::XDR xdr_file;
    std::vector<size_t> frame_indices;
    uint natoms_;
    GCoord box;
    double precision_;


  private:

    void scanFrames(void);
    bool readFrameHeader(Header& hdr);

    int sizeofint(int);
    int sizeofints(uint*, const uint);
    int decodebits(int*, uint);
    void decodeints(int*, const int, int, uint*, int*);


    int readCompressedCoords(std::vector<T>& fp)
    {
      int minint[3], maxint[3], *lip;
      int smallidx, minidx, maxidx;
      uint sizeint[3], sizesmall[3], bitsizeint[3], size3;
      int k, *buf1, *buf2, lsize, flag;
      int smallnum, smaller, larger, i, is_smaller, run;
      T precision, inv_precision;
      int tmp, *thiscoord,  prevcoord[3];
      unsigned int bitsize;
  
     
      if (!xdr_file.read(lsize))
        return(0);

      size3 = lsize * 3;

      /* Dont bother with compression for three atoms or less */
      if(lsize<=9) {
        float* tmp = new T[size3];
        xdr_file.read(tmp, size3);
        for (uint i=0; i<size3; ++i)
          fp.push_back(tmp[i]);
        delete[] tmp;
        return(lsize);
      }

      /* Compression-time if we got here. Read precision first */
      xdr_file.read(precision);
      precision_ = precision;
  
      int size3padded = static_cast<int>(size3 * 1.2);
      buf1 = new int[size3padded];
      buf2 = new int[size3padded];
      /* buf2[0-2] are special and do not contain actual data */
      buf2[0] = buf2[1] = buf2[2] = 0;
      xdr_file.read(minint, 3);
      xdr_file.read(maxint, 3);
  
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
	
      if (!xdr_file.read(smallidx)) {
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
  
      if (!xdr_file.read(buf2, 1)) {
        delete[] buf1;
        delete[] buf2;
        return(0);
      }

      if (!xdr_file.read(reinterpret_cast<char*>(&(buf2[3])), static_cast<uint>(buf2[0]))) {
        delete[] buf1;
        delete[] buf2;
        return(0);
      }
      buf2[0] = buf2[1] = buf2[2] = 0;
  
      inv_precision = 1.0 / precision;
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

              fp.push_back(prevcoord[0] * inv_precision);
              fp.push_back(prevcoord[1] * inv_precision);
              fp.push_back(prevcoord[2] * inv_precision);

            } else {
              prevcoord[0] = thiscoord[0];
              prevcoord[1] = thiscoord[1];
              prevcoord[2] = thiscoord[2];
            }
            fp.push_back(thiscoord[0] * inv_precision);
            fp.push_back(thiscoord[1] * inv_precision);
            fp.push_back(thiscoord[2] * inv_precision);

          }
        } else {
          fp.push_back(thiscoord[0] * inv_precision);
          fp.push_back(thiscoord[1] * inv_precision);
          fp.push_back(thiscoord[2] * inv_precision);
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
};


#endif
