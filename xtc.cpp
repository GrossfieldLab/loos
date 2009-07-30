#include <xtc.hpp>



namespace loos {


  template<typename T>
  int XTC<T>::sizeofint(int size) {
    int n = 0;
    while ( size > 0 ) {
      size >>= 1;
      ++n;
    }
    return(n);
  }

  
  template<typename T>
  int XTC<T>::sizeofints(uint* data, const uint n) {
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

  template<typename T>
  int XTC<T>::decodebits(int* buf, uint nbits) {

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

  template<typename T>
  void XTC<T>::decodeints(int* buf, const int num_of_ints, int num_of_bits,
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


  template<typename T>
  bool XTC<T>::readFrameHeader(XTC<T>::Header& hdr) {
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


  template<typename T>
  void XTC<T>::scanFrames(void) {
    frame_indices.clear();
    
    ifs()->seekg(0);

    Header h;

    while (! ifs()->eof()) {
      size_t pos = ifs()->tellg();
      bool ok = readFrameHeader(h);
      if (!ok)
        return;

      frame_indices.push_back(pos);
      if (natoms_ == 0)
        natoms_ = h.natoms;
      else if (natoms_ != h.natoms)
        throw(std::runtime_error("XTC frames have differing numbers of atoms"));

      uint block_size = xdr_file.block_size();
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


  template<typename T>
  void XTC<T>::seekFrameImpl(const uint i) {
    if (i >= frame_indices.size())
      throw(std::runtime_error("Trying to seek past the end of the file"));
    
    ifs()->clear();
    ifs()->seekg(frame_indices[i]);
  }



}
