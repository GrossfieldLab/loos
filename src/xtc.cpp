/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo, Alan Grossfield
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


/*
  NOTE:  This code is based on the xdrfile library authored by:
    David van der Spoel <spoel@gromacs.org>
    Erik Lindahl <lindahl@gromacs.org>
  and covered by the GLPL-v3 license
*/


#include <xtc.hpp>


namespace loos {


  // Globals...
  const int XTC::magicints[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
    80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024, 1290, // 31
    1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003, // 41
    16384, 20642, 26007, 32768, 41285, 52015, 65536,82570, 104031, // 50
    131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561, // 58
    832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021, // 65
    4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216 // 72
  };

  const int XTC::firstidx = 9;

  const int XTC::lastidx = 73;

  const int XTC::magic = 1995;

  const uint XTC::min_compressed_system_size = 9;
    


  // The following are largely from the xdrlib...
  int XTC::sizeofint(int size) {
    int n = 0;
    while ( size > 0 ) {
      size >>= 1;
      ++n;
    }
    return(n);
  }

  
  int XTC::sizeofints(uint* data, const uint n) {
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

  
  int XTC::decodebits(int* buf, uint nbits) {

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


  void XTC::decodeints(int* buf, const int nints, int nbits,
                       uint* sizes, int* nums) {
    int bytes[32];
    int i, j, num_of_bytes, p, num;
  
    bytes[1] = bytes[2] = bytes[3] = 0;
    num_of_bytes = 0;
    while (nbits > 8) {
      bytes[num_of_bytes++] = decodebits(buf, 8);
      nbits -= 8;
    }
    if (nbits > 0) {
      bytes[num_of_bytes++] = decodebits(buf, nbits);
    }
    for (i = nints-1; i > 0; i--) {
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



  // Coordinates are converted into GCoords and stored in the object's
  // coords_ vector

  bool XTC::readCompressedCoords(void)
  {
    int minint[3], maxint[3], *lip;
    int smallidx;
    uint sizeint[3], sizesmall[3], bitsizeint[3] = {0,0,0}, size3;
    int k, *buf1, *buf2, lsize, flag;
    int smallnum, smaller, i, is_smaller, run;
    xtc_t precision, inv_precision;
    int tmp, *thiscoord,  prevcoord[3];
    unsigned int bitsize;
  
     
    if (!xdr_file.read(lsize))
      return(false);

    size3 = lsize * 3;

    /* Dont bother with compression for three atoms or less */
    if(lsize<=9) {
      float* tmp = new xtc_t[size3];
      xdr_file.read(tmp, size3);
      for (uint i=0; i<size3; i += 3)
        coords_.push_back(GCoord(tmp[i], tmp[i+1], tmp[i+2]) * 10.0);
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
    tmp = smallidx-1;
    tmp = (firstidx>tmp) ? firstidx : tmp;
    smaller = magicints[tmp] / 2;
    smallnum = magicints[smallidx] / 2;
    sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;

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

            coords_.push_back(GCoord(prevcoord[0] * inv_precision,
                                prevcoord[1] * inv_precision,
                                prevcoord[2] * inv_precision) * 10.0);
          } else {
            prevcoord[0] = thiscoord[0];
            prevcoord[1] = thiscoord[1];
            prevcoord[2] = thiscoord[2];
          }
          coords_.push_back(GCoord(thiscoord[0] * inv_precision,
                              thiscoord[1] * inv_precision,
                              thiscoord[2] * inv_precision) * 10.0);
        }
      } else {
        coords_.push_back(GCoord(thiscoord[0] * inv_precision,
                            thiscoord[1] * inv_precision,
                            thiscoord[2] * inv_precision) * 10.0);
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



  bool XTC::readUncompressedCoords(void) 
  {
      uint lsize;
      
      if (!xdr_file.read(lsize))
	  return(false);
      
      uint size3 = lsize * 3;
      coords_ = std::vector<GCoord>(lsize);
      float* tmp_coords = new float[size3];
      uint n = xdr_file.read(tmp_coords, size3);
      if (n != size3)
	throw(FileReadError(_filename, "XTC Error: number of uncompressed coords read did not match number expected"));
      
      uint i = 0;
      for (uint j=0; j<lsize; ++j, i += 3)
	  coords_[j] = GCoord(tmp_coords[i], tmp_coords[i+1], tmp_coords[i+2]) * 10.0;
      
      delete[] tmp_coords;
      return(true);
  }
    
	      

  void XTC::updateGroupCoordsImpl(AtomicGroup& g) {

    for (AtomicGroup::iterator i = g.begin(); i != g.end(); ++i) {
      uint idx = (*i)->index();
      if (idx > natoms_)
        throw(TrajectoryError("updating group coords", _filename, "Atom index into trajectory frame is out of bounds"));
      (*i)->coords(coords_[idx]);
    }
    
    // XTC files *always* have a periodic box...
    g.periodicBox(box);
  }


  bool XTC::parseFrame(void) {
    if (ifs->eof())
      return(false);

    // First, clear out existing coords...  A read error after this
    // point will invalidate the current object's coord state

    coords_.clear();
    if (!readFrameHeader(current_header_))
      return(false);
    
    box = GCoord(current_header_.box[0], 
		 current_header_.box[4], 
		 current_header_.box[8]) * 10.0; // Convert to Angstroms
    if (natoms_ <= min_compressed_system_size)
	return(readUncompressedCoords());
    else
	return(readCompressedCoords());
  }


  bool XTC::readFrameHeader(XTC::Header& hdr) {
    int magic_no;
    int ok = xdr_file.read(magic_no);
    if (!ok)
      return(false);
    if (magic_no != magic) {
      std::ostringstream oss;
      oss << "Invalid XTC magic number (got " << magic_no << " but expected " << magic << ")";
      throw(FileReadError(_filename, oss.str()));
    }

    // Defer error-checks until the end...
    xdr_file.read(hdr.natoms);

    xdr_file.read(hdr.step);
    xdr_file.read(hdr.time);
    ok = xdr_file.read(hdr.box, 9);
    if (!ok)
      throw(FileReadError(_filename, "Problem reading XTC header"));

    return(true);
  }


  // Scan the trajectory file, skipping each compressed frame.  In the
  // process, we build up an index relating file-pos to frame index.
  // This permits fast seeking of indivual frames.
  void XTC::scanFrames(void) {
    frame_indices.clear();
    
    rewindImpl();

    Header h;
    
    while (! ifs->eof()) {
      size_t pos = ifs->tellg();

      bool ok = readFrameHeader(h);
      if (!ok) {
        rewindImpl();
        return;
      }

      frame_indices.push_back(pos);
      if (natoms_ == 0)
        natoms_ = h.natoms;
      else if (natoms_ != h.natoms)
        throw(FileOpenError(_filename, "XTC frames have differing numbers of atoms"));

      uint block_size = sizeof(internal::XDRReader::block_type);

      // Always update estimated timestep...
      if (h.step != 0)
	timestep_ = h.time / h.step;

      size_t offset = 0;
      uint nbytes = 0;

      if (natoms_ <= min_compressed_system_size) {
	  nbytes = natoms_ * 3 * sizeof(float);
	  uint dummy;
	  xdr_file.read(dummy);
	  if (dummy != natoms_)
	    throw(FileOpenError(_filename, "XTC small system vector size is not what was expected"));
      } else {
	  offset = 9 * block_size;
	  ifs->seekg(offset, std::ios_base::cur);
	  xdr_file.read(nbytes);
      }
      
      uint nblocks = nbytes / block_size;

      if (nbytes % block_size != 0)
        ++nblocks;   // round up
      offset = nblocks * block_size;
      ifs->seekg(offset, std::ios_base::cur);
    }

    // Catch-all for I/O errors
    if (ifs->fail())
      throw(FileOpenError(_filename, "Problem scanning XTC trajectory to build frame indices"));
    
  }


  void XTC::seekFrameImpl(const uint i) {
    if (i >= frame_indices.size())
      throw(FileError(_filename, "Requested XTC frame is out of range"));
    
    ifs->clear();
    ifs->seekg(frame_indices[i], std::ios_base::beg);
  }

}

