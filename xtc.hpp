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

    static const int magicints[];
    static const int firstidx, lastidx;

  public:
    explicit XTC(const std::string& s) : Trajectory(s), xdr_file(ifs()),natoms_(0) {
      scanFrames();
      coords_.reserve(natoms_);
      if (!parseFrame())
        throw(std::logic_error("Unable to read in the first frame"));
      cached_first = true;
    }
    explicit XTC(const char* p) : Trajectory(p), xdr_file(ifs()),natoms_(0) {
      scanFrames();
      coords_.reserve(natoms_);
      if (!parseFrame())
        throw(std::logic_error("Unable to read in the first frame"));
      cached_first = true;
    } 

    uint natoms(void) const { return(natoms_); }
    float timestep(void) const { return(0); }
    uint nframes(void) const { return(frame_indices.size()); }
    bool hasPeriodicBox(void) const { return(true); }
    GCoord periodicBox(void) const { return(box); }
    void updateGroupCoords(AtomicGroup& g) { }
    bool parseFrame(void) {
      if (ifs()->eof())
        return(false);

      coords_.clear();
      Header h;
      if (!readFrameHeader(h))
        return(false);

      box = GCoord(h.box[0], h.box[4], h.box[8]);
      return(readCompressedCoords(coords_));
    }

    std::vector<GCoord> coords(void) {
      std::vector<GCoord> result;
      for (uint i=0; i<coords_.size(); i += 3)
        result.push_back(GCoord(coords_[i], coords_[i+1], coords_[i+2]));

      return(result);
    }
      

    void seekNextFrameImpl(void) { }
    void rewindImpl(void) { ifs()->clear(); ifs()->seekg(0); }

  private:

    internal::XDR xdr_file;
    std::vector<size_t> frame_indices;
    uint natoms_;
    GCoord box;
    double precision_;
    std::vector<T> coords_;


  private:


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


  

    bool readFrameHeader(XTC<T>::Header& hdr) {
      int magic;
      int ok = xdr_file.read(magic);
      if (!ok)
        return(false);
      if (magic != 1995)
        throw(std::runtime_error("Invalid XTC magic number"));

      xdr_file.read(hdr.natoms);

      xdr_file.read(hdr.step);
      xdr_file.read(hdr.time);
      ok = xdr_file.read(hdr.box, 9);
      if (!ok)
        throw(std::runtime_error("Error while reading XTC frame header"));

      return(true);
    }


    void scanFrames(void) {
      frame_indices.clear();
    
      rewindImpl();

      Header h;

      while (! ifs()->eof()) {
        size_t pos = ifs()->tellg();
        bool ok = readFrameHeader(h);
        if (!ok) {
          rewindImpl();
          return;
        }

        frame_indices.push_back(pos);
        if (natoms_ == 0)
          natoms_ = h.natoms;
        else if (natoms_ != h.natoms)
          throw(std::runtime_error("XTC frames have differing numbers of atoms"));

        uint block_size = sizeof(internal::XDR::block_type);
        size_t offset = 9 * block_size;
        ifs()->seekg(offset, std::ios_base::cur);
        uint nbytes;
        xdr_file.read(nbytes);
        uint nblocks = nbytes / block_size;
        if (nbytes % block_size != 0)
          ++nblocks;   // round up
        offset = nblocks * block_size;
        ifs()->seekg(offset, std::ios_base::cur);
      }
    }


    void seekFrameImpl(const uint i) {
      if (i >= frame_indices.size())
        throw(std::runtime_error("Trying to seek past the end of the file"));
    
      ifs()->clear();
      ifs()->seekg(frame_indices[i], std::ios_base::beg);
    }




    bool readCompressedCoords(std::vector<T>& fp)
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
        return(false);

      size3 = lsize * 3;

      /* Dont bother with compression for three atoms or less */
      if(lsize<=9) {
        float* tmp = new T[size3];
        xdr_file.read(tmp, size3);
        for (uint i=0; i<size3; ++i)
          fp.push_back(tmp[i]);
        delete[] tmp;
        return(true);
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
        return(false);
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
        return(false);
      }

      if (!xdr_file.read(reinterpret_cast<char*>(&(buf2[3])), static_cast<uint>(buf2[0]))) {
        delete[] buf1;
        delete[] buf2;
        return(false);
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
      return(true);
    }


  };


  template<typename T>
  const int XTC<T>::magicints[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
    80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024, 1290, // 31
    1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003, // 41
    16384, 20642, 26007, 32768, 41285, 52015, 65536,82570, 104031, // 50
    131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561, // 58
    832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021, // 65
    4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216 // 72
  };

  template<typename T>
  const int XTC<T>::firstidx = 9;

  template<typename T>
  const int XTC<T>::lastidx = 73;


};


#endif
