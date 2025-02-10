/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2014, Tod D. Romo, Alan Grossfield
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


#include <xtcwriter.hpp>
#include <xtc.hpp>

namespace loos 
{
  

  // -- The following is adapted from xdrfile-1.1b --
  // Copyright (c) Erik Lindahl, David van der Spoel 2003,2004.
  // Coordinate compression (c) by Frans van Hoesel. 

  // http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library

  /******************************************************************
                   GNU LESSER GENERAL PUBLIC LICENSE
                       Version 3, 29 June 2007

 Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
 Everyone is permitted to copy and distribute verbatim copies
 of this license document, but changing it is not allowed.


  This version of the GNU Lesser General Public License incorporates
the terms and conditions of version 3 of the GNU General Public
License, supplemented by the additional permissions listed below.

  0. Additional Definitions.

  As used herein, "this License" refers to version 3 of the GNU Lesser
General Public License, and the "GNU GPL" refers to version 3 of the GNU
General Public License.

  "The Library" refers to a covered work governed by this License,
other than an Application or a Combined Work as defined below.

  An "Application" is any work that makes use of an interface provided
by the Library, but which is not otherwise based on the Library.
Defining a subclass of a class defined by the Library is deemed a mode
of using an interface provided by the Library.

  A "Combined Work" is a work produced by combining or linking an
Application with the Library.  The particular version of the Library
with which the Combined Work was made is also called the "Linked
Version".

  The "Minimal Corresponding Source" for a Combined Work means the
Corresponding Source for the Combined Work, excluding any source code
for portions of the Combined Work that, considered in isolation, are
based on the Application, and not on the Linked Version.

  The "Corresponding Application Code" for a Combined Work means the
object code and/or source code for the Application, including any data
and utility programs needed for reproducing the Combined Work from the
Application, but excluding the System Libraries of the Combined Work.

  1. Exception to Section 3 of the GNU GPL.

  You may convey a covered work under sections 3 and 4 of this License
without being bound by section 3 of the GNU GPL.

  2. Conveying Modified Versions.

  If you modify a copy of the Library, and, in your modifications, a
facility refers to a function or data to be supplied by an Application
that uses the facility (other than as an argument passed when the
facility is invoked), then you may convey a copy of the modified
version:

   a) under this License, provided that you make a good faith effort to
   ensure that, in the event an Application does not supply the
   function or data, the facility still operates, and performs
   whatever part of its purpose remains meaningful, or

   b) under the GNU GPL, with none of the additional permissions of
   this License applicable to that copy.

  3. Object Code Incorporating Material from Library Header Files.

  The object code form of an Application may incorporate material from
a header file that is part of the Library.  You may convey such object
code under terms of your choice, provided that, if the incorporated
material is not limited to numerical parameters, data structure
layouts and accessors, or small macros, inline functions and templates
(ten or fewer lines in length), you do both of the following:

   a) Give prominent notice with each copy of the object code that the
   Library is used in it and that the Library and its use are
   covered by this License.

   b) Accompany the object code with a copy of the GNU GPL and this license
   document.

  4. Combined Works.

  You may convey a Combined Work under terms of your choice that,
taken together, effectively do not restrict modification of the
portions of the Library contained in the Combined Work and reverse
engineering for debugging such modifications, if you also do each of
the following:

   a) Give prominent notice with each copy of the Combined Work that
   the Library is used in it and that the Library and its use are
   covered by this License.

   b) Accompany the Combined Work with a copy of the GNU GPL and this license
   document.

   c) For a Combined Work that displays copyright notices during
   execution, include the copyright notice for the Library among
   these notices, as well as a reference directing the user to the
   copies of the GNU GPL and this license document.

   d) Do one of the following:

       0) Convey the Minimal Corresponding Source under the terms of this
       License, and the Corresponding Application Code in a form
       suitable for, and under terms that permit, the user to
       recombine or relink the Application with a modified version of
       the Linked Version to produce a modified Combined Work, in the
       manner specified by section 6 of the GNU GPL for conveying
       Corresponding Source.

       1) Use a suitable shared library mechanism for linking with the
       Library.  A suitable mechanism is one that (a) uses at run time
       a copy of the Library already present on the user's computer
       system, and (b) will operate properly with a modified version
       of the Library that is interface-compatible with the Linked
       Version.

   e) Provide Installation Information, but only if you would otherwise
   be required to provide such information under section 6 of the
   GNU GPL, and only to the extent that such information is
   necessary to install and execute a modified version of the
   Combined Work produced by recombining or relinking the
   Application with a modified version of the Linked Version. (If
   you use option 4d0, the Installation Information must accompany
   the Minimal Corresponding Source and Corresponding Application
   Code. If you use option 4d1, you must provide the Installation
   Information in the manner specified by section 6 of the GNU GPL
   for conveying Corresponding Source.)

  5. Combined Libraries.

  You may place library facilities that are a work based on the
Library side by side in a single library together with other library
facilities that are not Applications and are not covered by this
License, and convey such a combined library under terms of your
choice, if you do both of the following:

   a) Accompany the combined library with a copy of the same work based
   on the Library, uncombined with any other library facilities,
   conveyed under the terms of this License.

   b) Give prominent notice with the combined library that part of it
   is a work based on the Library, and explaining where to find the
   accompanying uncombined form of the same work.

  6. Revised Versions of the GNU Lesser General Public License.

  The Free Software Foundation may publish revised and/or new versions
of the GNU Lesser General Public License from time to time. Such new
versions will be similar in spirit to the present version, but may
differ in detail to address new problems or concerns.

  Each version is given a distinguishing version number. If the
Library as you received it specifies that a certain numbered version
of the GNU Lesser General Public License "or any later version"
applies to it, you have the option of following the terms and
conditions either of that published version or of any later version
published by the Free Software Foundation. If the Library as you
received it does not specify a version number of the GNU Lesser
General Public License, you may choose any version of the GNU Lesser
General Public License ever published by the Free Software Foundation.

  If the Library as you received it specifies that a proxy can decide
whether future versions of the GNU Lesser General Public License shall
apply, that proxy's public statement of acceptance of any version is
permanent authorization for you to choose that version for the
Library.

******************************************************************/

  // TDR - This needs to move?
  const int XTCWriter::magicints[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
    80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024, 1290,
    1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003, 
    16384, 20642, 26007, 32768, 41285, 52015, 65536,82570, 104031, 
    131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561, 
    832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021, 
    4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216 
  };


  const int XTCWriter::firstidx = 9;
  const int XTCWriter::lastidx = sizeof(magicints) / sizeof(*magicints);

  const int XTCWriter::DIM = 3;

  /* Internal support routines for reading/writing compressed coordinates 
   * sizeofint - calculate smallest number of bits necessary
   * to represent a certain integer.
   */
  int XTCWriter::sizeofint(const int size) const {
    int num = 1;
    int num_of_bits = 0;
    
    while (size >= num && num_of_bits < 32) {
      num_of_bits++;
      num <<= 1;
    }
    return num_of_bits;
  }


  /*
   * sizeofints - calculate 'bitsize' of compressed ints
   *
   * given a number of small unsigned integers and the maximum value
   * return the number of bits needed to read or write them with the
   * routines encodeints/decodeints. You need this parameter when
   * calling those routines. 
   * (However, in some cases we can just use the variable 'smallidx' 
   * which is the exact number of bits, and them we dont need to call
   * this routine).
   */
  int XTCWriter::sizeofints(const int num_of_ints, const unsigned int sizes[]) const
  {
    int i, num;
    unsigned int num_of_bytes, num_of_bits, bytes[32], bytecnt, tmp;
    num_of_bytes = 1;
    bytes[0] = 1;
    num_of_bits = 0;
    for (i=0; i < num_of_ints; i++) {   
      tmp = 0;
      for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++) {
        tmp = bytes[bytecnt] * sizes[i] + tmp;
        bytes[bytecnt] = tmp & 0xff;
        tmp >>= 8;
      }
      while (tmp != 0) {
        bytes[bytecnt++] = tmp & 0xff;
        tmp >>= 8;
      }
      num_of_bytes = bytecnt;
    }
    num = 1;
    num_of_bytes--;
    while (bytes[num_of_bytes] >= num) {
      num_of_bits++;
      num *= 2;
    }
    return(num_of_bits + num_of_bytes * 8);
  }
    

  /*
   * encodebits - encode num into buf using the specified number of bits
   *
   * This routines appends the value of num to the bits already present in
   * the array buf. You need to give it the number of bits to use and you had
   * better make sure that this number of bits is enough to hold the value.
   * Num must also be positive.
   */
  void XTCWriter::encodebits(int buf[], int num_of_bits, const int num) const
  {
    
    unsigned int cnt, lastbyte;
    int lastbits;
    unsigned char * cbuf;
    
    cbuf = ((unsigned char *)buf) + 3 * sizeof(*buf);
    cnt = (unsigned int) buf[0];
    lastbits = buf[1];
    lastbyte =(unsigned int) buf[2];
    while (num_of_bits >= 8) {
      lastbyte = (lastbyte << 8) | ((num >> (num_of_bits -8)) /* & 0xff*/);
      cbuf[cnt++] = lastbyte >> lastbits;
      num_of_bits -= 8;
    }
    if (num_of_bits > 0) {
      lastbyte = (lastbyte << num_of_bits) | num;
      lastbits += num_of_bits;
      if (lastbits >= 8) {
        lastbits -= 8;
        cbuf[cnt++] = lastbyte >> lastbits;
      }
    }
    buf[0] = cnt;
    buf[1] = lastbits;
    buf[2] = lastbyte;
    if (lastbits>0) {
      cbuf[cnt] = lastbyte << (8 - lastbits);
    }
  }

  /*
   * encodeints - encode a small set of small integers in compressed format
   *
   * this routine is used internally by xdr3dfcoord, to encode a set of
   * small integers to the buffer for writing to a file.
   * Multiplication with fixed (specified maximum) sizes is used to get
   * to one big, multibyte integer. Allthough the routine could be
   * modified to handle sizes bigger than 16777216, or more than just
   * a few integers, this is not done because the gain in compression
   * isn't worth the effort. Note that overflowing the multiplication
   * or the byte buffer (32 bytes) is unchecked and whould cause bad results.
   * THese things are checked in the calling routines, so make sure not
   * to remove those checks...
   */
 
  void XTCWriter::encodeints(int buf[], const int num_of_ints, const int num_of_bits,
                             const unsigned int sizes[], const unsigned int nums[]) const
  {

    int i;
    unsigned int bytes[32], num_of_bytes, bytecnt, tmp;

    tmp = nums[0];
    num_of_bytes = 0;
    do {
      bytes[num_of_bytes++] = tmp & 0xff;
      tmp >>= 8;
    } while (tmp != 0);

    for (i = 1; i < num_of_ints; i++) {
      if (nums[i] >= sizes[i]) {
        std::ostringstream oss;
        oss << boost::format("Major breakdown in XTCWriter::encodeints() - num %u doesn't match size %u")
          % nums[i]
          % sizes[i];
        throw(LOOSError(oss.str()));
      }
      /* use one step multiply */    
      tmp = nums[i];
      for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++) {
        tmp = bytes[bytecnt] * sizes[i] + tmp;
        bytes[bytecnt] = tmp & 0xff;
        tmp >>= 8;
      }
      while (tmp != 0) {
        bytes[bytecnt++] = tmp & 0xff;
        tmp >>= 8;
      }
      num_of_bytes = bytecnt;
    }
    if (num_of_bits >= num_of_bytes * 8) {
      for (i = 0; i < num_of_bytes; i++) {
        encodebits(buf, 8, bytes[i]);
      }
      encodebits(buf, num_of_bits - num_of_bytes * 8, 0);
    } 
    else {
      for (i = 0; i < num_of_bytes-1; i++){
        encodebits(buf, 8, bytes[i]);
      }
      encodebits(buf, num_of_bits- (num_of_bytes -1) * 8, bytes[i]);
    }
  }





  void XTCWriter::writeCompressedCoordsFloat(float* ptr, int size, float precision) 
  {
    int minint[3], maxint[3], mindiff, *lip, diff;
    int lint1, lint2, lint3, oldlint1, oldlint2, oldlint3, smallidx;
    int minidx, maxidx;
    unsigned sizeint[3], sizesmall[3], bitsizeint[3], *luip;
    int k;
    int smallnum, smaller, larger, i, j, is_small, is_smaller, run, prevrun;
    float *lfp, lf;
    int tmp, tmpsum, *thiscoord,  prevcoord[3];
    unsigned int tmpcoord[30];
    unsigned int bitsize;

    uint size3 = size * 3;
  
    bitsizeint[0] = 0;
    bitsizeint[1] = 0;
    bitsizeint[2] = 0;

    allocateBuffers(size);
    if (!xdr.write(size))
      throw(FileWriteError(_filename, "Could not write size to XTC file"));

    /* Dont bother with compression for three atoms or less */
    if(size<=9) 
    {
      xdr.write(ptr, size3);
      return;
    }
    /* Compression-time if we got here. Write precision first */
    if (precision <= 0)
      precision = 1000;

    xdr.write(precision);
    /* buf2[0-2] are special and do not contain actual data */
    buf2[0] = buf2[1] = buf2[2] = 0;
    minint[0] = minint[1] = minint[2] = INT_MAX;
    maxint[0] = maxint[1] = maxint[2] = INT_MIN;
    prevrun = -1;
    lfp = ptr;
    lip = buf1;
    mindiff = INT_MAX;
    oldlint1 = oldlint2 = oldlint3 = 0;
    while(lfp < ptr + size3 )
    {
      /* find nearest integer */
      if (*lfp >= 0.0)
        lf = *lfp * precision + 0.5;
      else
        lf = *lfp * precision - 0.5;
      if (fabs(lf) > INT_MAX-2) 
      {
        /* scaling would cause overflow */
        throw(LOOSError("Internal overflow compressing coordinates...check input model coordinates (#1)"));
      }
      lint1 = lf;
      if (lint1 < minint[0]) minint[0] = lint1;
      if (lint1 > maxint[0]) maxint[0] = lint1;
      *lip++ = lint1;
      lfp++;
      if (*lfp >= 0.0)
        lf = *lfp * precision + 0.5;
      else
        lf = *lfp * precision - 0.5;
      if (fabs(lf) > INT_MAX-2)
      {
        /* scaling would cause overflow */
        throw(LOOSError("Internal overflow compressing coordinates...check input model coordinates (#2)"));
      }
      lint2 = lf;
      if (lint2 < minint[1]) minint[1] = lint2;
      if (lint2 > maxint[1]) maxint[1] = lint2;
      *lip++ = lint2;
      lfp++;
      if (*lfp >= 0.0)
        lf = *lfp * precision + 0.5;
      else
        lf = *lfp * precision - 0.5;

      // *** TDR - This is not actually used
      // if (fabs(lf) > INT_MAX-2) 
      // {
      //        errval=0;      
      // }
      lint3 = lf;
      if (lint3 < minint[2]) minint[2] = lint3;
      if (lint3 > maxint[2]) maxint[2] = lint3;
      *lip++ = lint3;
      lfp++;
      diff = abs(oldlint1-lint1)+abs(oldlint2-lint2)+abs(oldlint3-lint3);
      if (diff < mindiff && lfp > ptr + 3)
        mindiff = diff;
      oldlint1 = lint1;
      oldlint2 = lint2;
      oldlint3 = lint3;
    }  
    xdr.write(minint, 3);
    xdr.write(maxint, 3);
  
    if ((float)maxint[0] - (float)minint[0] >= INT_MAX-2 ||
        (float)maxint[1] - (float)minint[1] >= INT_MAX-2 ||
        (float)maxint[2] - (float)minint[2] >= INT_MAX-2) {
      /* turning value to unsigned by subtracting minint
       * would cause overflow
       */
      throw(LOOSError("Internal overflow compressing internal coordinates...check input model coordinates (#3)"));
    }
    sizeint[0] = maxint[0] - minint[0]+1;
    sizeint[1] = maxint[1] - minint[1]+1;
    sizeint[2] = maxint[2] - minint[2]+1;
  
    /* check if one of the sizes is to big to be multiplied */
    if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff)
    {
      bitsizeint[0] = sizeofint(sizeint[0]);
      bitsizeint[1] = sizeofint(sizeint[1]);
      bitsizeint[2] = sizeofint(sizeint[2]);
      bitsize = 0; /* flag the use of large sizes */
    }
    else
    {
      bitsize = sizeofints(3, sizeint);
    }
    lip = buf1;
    luip = (unsigned int *) buf1;
    smallidx = firstidx;
    while (smallidx < lastidx && magicints[smallidx] < mindiff)
    {
      smallidx++;
    }
    xdr.write(smallidx);
    tmp=smallidx+8;
    maxidx = (lastidx<tmp) ? lastidx : tmp;
    minidx = maxidx - 8; /* often this equal smallidx */
    tmp=smallidx-1;
    tmp= (firstidx>tmp) ? firstidx : tmp;
    smaller = magicints[tmp] / 2;
    smallnum = magicints[smallidx] / 2;
    sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
    larger = magicints[maxidx] / 2;
    i = 0;
    while (i < size) 
    {
      is_small = 0;
      thiscoord = (int *)(luip) + i * 3;
      if (smallidx < maxidx && i >= 1 &&
          abs(thiscoord[0] - prevcoord[0]) < larger &&
          abs(thiscoord[1] - prevcoord[1]) < larger &&
          abs(thiscoord[2] - prevcoord[2]) < larger) {
        is_smaller = 1;
      } 
      else if (smallidx > minidx) 
      {
        is_smaller = -1;
      }
      else
      {
        is_smaller = 0;
      }
      if (i + 1 < size) 
      {
        if (abs(thiscoord[0] - thiscoord[3]) < smallnum &&
            abs(thiscoord[1] - thiscoord[4]) < smallnum &&
            abs(thiscoord[2] - thiscoord[5]) < smallnum) 
        {
          /* interchange first with second atom for better
           * compression of water molecules
           */
          tmp = thiscoord[0]; thiscoord[0] = thiscoord[3];
          thiscoord[3] = tmp;
          tmp = thiscoord[1]; thiscoord[1] = thiscoord[4];
          thiscoord[4] = tmp;
          tmp = thiscoord[2]; thiscoord[2] = thiscoord[5];
          thiscoord[5] = tmp;
          is_small = 1;
        } 
      }
      tmpcoord[0] = thiscoord[0] - minint[0];
      tmpcoord[1] = thiscoord[1] - minint[1];
      tmpcoord[2] = thiscoord[2] - minint[2];
      if (bitsize == 0) 
      {
        encodebits(buf2, bitsizeint[0], tmpcoord[0]);
        encodebits(buf2, bitsizeint[1], tmpcoord[1]);
        encodebits(buf2, bitsizeint[2], tmpcoord[2]);
      } 
      else
      {
        encodeints(buf2, 3, bitsize, sizeint, tmpcoord);
      }
      prevcoord[0] = thiscoord[0];
      prevcoord[1] = thiscoord[1];
      prevcoord[2] = thiscoord[2];
      thiscoord = thiscoord + 3;
      i++;

      run = 0;
      if (is_small == 0 && is_smaller == -1)
        is_smaller = 0;
      while (is_small && run < 8*3)
      {
        tmpsum=0;
        for(j=0;j<3;j++) 
        {
          tmp=thiscoord[j] - prevcoord[j];
          tmpsum+=tmp*tmp;
        }
        if (is_smaller == -1 && tmpsum >= smaller * smaller)
        {
          is_smaller = 0;
        }
      
        tmpcoord[run++] = thiscoord[0] - prevcoord[0] + smallnum;
        tmpcoord[run++] = thiscoord[1] - prevcoord[1] + smallnum;
        tmpcoord[run++] = thiscoord[2] - prevcoord[2] + smallnum;
      
        prevcoord[0] = thiscoord[0];
        prevcoord[1] = thiscoord[1];
        prevcoord[2] = thiscoord[2];
      
        i++;
        thiscoord = thiscoord + 3;
        is_small = 0;
        if (i < size &&
            abs(thiscoord[0] - prevcoord[0]) < smallnum &&
            abs(thiscoord[1] - prevcoord[1]) < smallnum &&
            abs(thiscoord[2] - prevcoord[2]) < smallnum)
        {
          is_small = 1;
        }
      }
      if (run != prevrun || is_smaller != 0) 
      {
        prevrun = run;
        encodebits(buf2, 1, 1); /* flag the change in run-length */
        encodebits(buf2, 5, run+is_smaller+1);
      } 
      else 
      {
        encodebits(buf2, 1, 0); /* flag the fact that runlength did not change */
      }
      for (k=0; k < run; k+=3) 
      {
        encodeints(buf2, 3, smallidx, sizesmall, &tmpcoord[k]); 
      }
      if (is_smaller != 0) 
      {
        smallidx += is_smaller;
        if (is_smaller < 0) 
        {
          smallnum = smaller;
          smaller = magicints[smallidx-1] / 2;
        } 
        else 
        {
          smaller = smallnum;
          smallnum = magicints[smallidx] / 2;
        }
        sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
      }   
    }
    if (buf2[1] != 0) buf2[0]++;
    xdr.write(buf2[0]);
    tmp=xdr.write((char *)&(buf2[3]),(unsigned int)buf2[0]);
    if(tmp!=(unsigned int)buf2[0])
      throw(FileWriteError(_filename, "Error while writing compressed coordinates to XTC file"));
  }



  // -- End of code from xdrfile library --


  // Handle allocation of buffers (would be handle by system xdr lib)
  void XTCWriter::allocateBuffers(const size_t size) {
    size_t size3 = size * 3;
    if (size3 > buf1size) {
      if (buf1)
        delete[] buf1;
      if (buf2)
        delete[] buf2;

      buf1 = new int[size3];
      size_t size3plus = size3 * 1.2;
      buf2 = new int[size3plus];

      buf1size = size3;
      buf2size = size3 * 1.2;
    }
  }


  // Write a frame header
  void XTCWriter::writeHeader(const int natoms, const int step, const float time) {
    int magic = 1995;

    xdr.write(magic);
    xdr.write(natoms);
    xdr.write(step);
    xdr.write(time);
  }


  // Write a periodic box, translating from A to nm
  void XTCWriter::writeBox(const GCoord& box) {
    float outbox[DIM*DIM];
    for (uint i=0; i < DIM*DIM; ++i)
      outbox[i] = 0.0;
    outbox[0] = box[0] / 10.0;
    outbox[4] = box[1] / 10.0;
    outbox[8] = box[2] / 10.0; 

    xdr.write(outbox, DIM*DIM);
  }

  

  // Write a frame, converting units from A to nm.  Will allocate a temp array to hold coords...
  void XTCWriter::writeFrame(const AtomicGroup& model, const uint step, const double time) {

    writeHeader(model.size(), step, time);
    writeBox(model.periodicBox());
    uint n = model.size();

    if (n > crds_size_) {
      delete[] crds_;
      crds_ = new float[n * 3];
      crds_size_ = n;
    }

    for (uint i=0,k=0; i<n; ++i) {
      GCoord c = model[i]->coords();
      crds_[k++] = c.x() / 10.0;       // Convert to nm
      crds_[k++] = c.y() / 10.0;
      crds_[k++] = c.z() / 10.0;
    }
    writeCompressedCoordsFloat(crds_, n, precision_);

    ++current_;
  }



  void XTCWriter::writeFrame(const AtomicGroup& model) {
    writeFrame(model, step_, dt_ * step_);
    step_ += steps_per_frame_;
  }


  // Read existing XTC to get frame count...
  void XTCWriter::prepareToAppend() {
    stream_->seekg(0);
    XTC xtc(*stream_);
    current_ = xtc.nframes();
    stream_->seekp(0, std::ios_base::end);
  }

};
