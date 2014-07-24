// Nascent xtc writer class

#include <xtcwriter.hpp>

namespace loos 
{
  

  // TDR - This needs to move
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


  /* Internal support routines for reading/writing compressed coordinates 
   * sizeofint - calculate smallest number of bits necessary
   * to represent a certain integer.
   */
  int XTCWriter::sizeofint(const int size) const {
    unsigned int num = 1;
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
	throw(std::runtime_error(oss.str()));
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





  int XTCWriter::writeCompressedCoordsFloat(float* ptr, int size, float precision) 
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
    int errval=1;
    unsigned int bitsize;

    uint size3 = size * 3;
  
    bitsizeint[0] = 0;
    bitsizeint[1] = 0;
    bitsizeint[2] = 0;

    allocateBuffers(size);
    if (!xdr.write(&size))
      return -1; /* return if we could not write size */
    /* Dont bother with compression for three atoms or less */
    if(size<=9) 
    {
      return(xdr.write(ptr, size3));
      /* return number of coords, not floats */
    }
    /* Compression-time if we got here. Write precision first */
    if (precision <= 0)
      precision = 1000;

    xdr.write(&precision);
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
	fprintf(stderr,"Internal overflow compressing coordinates.\n");
	errval=0;
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
	fprintf(stderr,"Internal overflow compressing coordinates.\n");
	errval=0;
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
      if (fabs(lf) > INT_MAX-2) 
      {
	errval=0;      
      }
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
      /* turning value in unsigned by subtracting minint
       * would cause overflow
       */
      fprintf(stderr,"Internal overflow compressing coordinates.\n");
      errval=0;
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
    xdr.write(&smallidx);
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
    xdr.write(buf2);
    tmp=xdr.write((char *)&(buf2[3]),(unsigned int)buf2[0]);
    if(tmp==(unsigned int)buf2[0])
      return size;
    else
      return -1;
  }




#if defined(FOOBARFLUBBER)

  int XTCWriter::writeCompressedCoordsDouble(double   *ptr,
					     int      size,
					     double    precision)
  {
    int minint[3], maxint[3], mindiff, *lip, diff;
    int lint1, lint2, lint3, oldlint1, oldlint2, oldlint3, smallidx;
    int minidx, maxidx;
    unsigned sizeint[3], sizesmall[3], bitsizeint[3], size3, *luip;
    int k, *buf1, *buf2;
    int smallnum, smaller, larger, i, j, is_small, is_smaller, run, prevrun;
    double *lfp;
    float float_prec, lf,tmpdata[30];
    int tmp, tmpsum, *thiscoord,  prevcoord[3];
    unsigned int tmpcoord[30];
    int errval=1;
    unsigned int bitsize;
  
    bitsizeint[0] = 0;
    bitsizeint[1] = 0;
    bitsizeint[2] = 0;

    if(xfp==NULL)
      return -1;
    size3=3*size;
    if(size3>xfp->buf1size) {
      if((xfp->buf1=(int *)malloc(sizeof(int)*size3))==NULL) {
	fprintf(stderr,"Cannot allocate memory for compressing coordinates.\n");
	return -1;
      }
      xfp->buf1size=size3;
      xfp->buf2size=size3*1.2;
      if((xfp->buf2=(int *)malloc(sizeof(int)*xfp->buf2size))==NULL) {
	fprintf(stderr,"Cannot allocate memory for compressing coordinates.\n");
	return -1;
      }
    }
    if(xdrfile_write_int(&size,1,xfp)==0)
      return -1; /* return if we could not write size */
    /* Dont bother with compression for three atoms or less */
    // TDR - this is one diff with float...copy doubles into float and then write out...
    if(size<=9) {
      for(i=0;i<9*3;i++)
	tmpdata[i]=ptr[i];
      return xdrfile_write_float(tmpdata,size3,xfp)/3;
      /* return number of coords, not floats */
    }
    /* Compression-time if we got here. Write precision first */
    if (precision <= 0)
      precision = 1000;

    // TDR - another diff with float...precision is always written as floating
    float_prec=precision;
    xdrfile_write_float(&float_prec,1,xfp);
    /* avoid repeated pointer dereferencing. */
    buf1=xfp->buf1; 
    buf2=xfp->buf2;
    /* buf2[0-2] are special and do not contain actual data */
    buf2[0] = buf2[1] = buf2[2] = 0;
    minint[0] = minint[1] = minint[2] = INT_MAX;
    maxint[0] = maxint[1] = maxint[2] = INT_MIN;
    prevrun = -1;
    lfp = ptr;
    lip = buf1;
    mindiff = INT_MAX;
    oldlint1 = oldlint2 = oldlint3 = 0;
    while(lfp < ptr + size3 ) {
      /* find nearest integer */
      if (*lfp >= 0.0)
	lf = (float)*lfp * float_prec + 0.5;
      else
	lf = (float)*lfp * float_prec - 0.5;
      if (fabs(lf) > INT_MAX-2) {
	/* scaling would cause overflow */
	fprintf(stderr,"Internal overflow compressing coordinates.\n");
	errval=0;
      }
      lint1 = lf;
      if (lint1 < minint[0]) minint[0] = lint1;
      if (lint1 > maxint[0]) maxint[0] = lint1;
      *lip++ = lint1;
      lfp++;
      if (*lfp >= 0.0)
	lf = (float)*lfp * float_prec + 0.5;
      else
	lf = (float)*lfp * float_prec - 0.5;
      if (fabs(lf) > INT_MAX-2) {
	/* scaling would cause overflow */
	fprintf(stderr,"Internal overflow compressing coordinates.\n");
	errval=0;
      }
      lint2 = lf;
      if (lint2 < minint[1]) minint[1] = lint2;
      if (lint2 > maxint[1]) maxint[1] = lint2;
      *lip++ = lint2;
      lfp++;
      if (*lfp >= 0.0)
	lf = (float)*lfp * float_prec + 0.5;
      else
	lf = (float)*lfp * float_prec - 0.5;
      if (fabs(lf) > INT_MAX-2) {
	errval=0;      
      }
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
    xdrfile_write_int(minint,3,xfp);
    xdrfile_write_int(maxint,3,xfp);
  
    if ((float)maxint[0] - (float)minint[0] >= INT_MAX-2 ||
	(float)maxint[1] - (float)minint[1] >= INT_MAX-2 ||
	(float)maxint[2] - (float)minint[2] >= INT_MAX-2) {
      /* turning value in unsigned by subtracting minint
       * would cause overflow
       */
      fprintf(stderr,"Internal overflow compressing coordinates.\n");
      errval=0;
    }
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
      bitsize = sizeofints(3, sizeint);
    }
    lip = buf1;
    luip = (unsigned int *) buf1;
    smallidx = firstidx;
    while (smallidx < lastidx && magicints[smallidx] < mindiff) {
      smallidx++;
    }
    xdrfile_write_int(&smallidx,1,xfp);
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
    while (i < size) {
      is_small = 0;
      thiscoord = (int *)(luip) + i * 3;
      if (smallidx < maxidx && i >= 1 &&
	  abs(thiscoord[0] - prevcoord[0]) < larger &&
	  abs(thiscoord[1] - prevcoord[1]) < larger &&
	  abs(thiscoord[2] - prevcoord[2]) < larger) {
	is_smaller = 1;
      } else if (smallidx > minidx) {
	is_smaller = -1;
      } else {
	is_smaller = 0;
      }
      if (i + 1 < size) {
	if (abs(thiscoord[0] - thiscoord[3]) < smallnum &&
	    abs(thiscoord[1] - thiscoord[4]) < smallnum &&
	    abs(thiscoord[2] - thiscoord[5]) < smallnum) {
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
      if (bitsize == 0) {
	encodebits(buf2, bitsizeint[0], tmpcoord[0]);
	encodebits(buf2, bitsizeint[1], tmpcoord[1]);
	encodebits(buf2, bitsizeint[2], tmpcoord[2]);
      } else {
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
      while (is_small && run < 8*3) {
	tmpsum=0;
	for(j=0;j<3;j++) {
	  tmp=thiscoord[j] - prevcoord[j];
	  tmpsum+=tmp*tmp;
	}
	if (is_smaller == -1 && tmpsum >= smaller * smaller) {
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
	    abs(thiscoord[2] - prevcoord[2]) < smallnum) {
	  is_small = 1;
	}
      }
      if (run != prevrun || is_smaller != 0) {
	prevrun = run;
	encodebits(buf2, 1, 1); /* flag the change in run-length */
	encodebits(buf2, 5, run+is_smaller+1);
      } else {
	encodebits(buf2, 1, 0); /* flag the fact that runlength did not change */
      }
      for (k=0; k < run; k+=3) {
	encodeints(buf2, 3, smallidx, sizesmall, &tmpcoord[k]);	
      }
      if (is_smaller != 0) {
	smallidx += is_smaller;
	if (is_smaller < 0) {
	  smallnum = smaller;
	  smaller = magicints[smallidx-1] / 2;
	} else {
	  smaller = smallnum;
	  smallnum = magicints[smallidx] / 2;
	}
	sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
      }   
    }
    if (buf2[1] != 0) buf2[0]++;
    xdrfile_write_int(buf2,1,xfp); /* buf2[0] holds the length in bytes */
    tmp=xdrfile_write_opaque((char *)&(buf2[3]),(unsigned int)buf2[0],xfp);
    if(tmp==(unsigned int)buf2[0])
      return size;
    else
      return -1; 
  }

#endif


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



  void writeHeader(const int natoms, const int step, const float time) {
    int magic = 1995;

    xdr.write(&magic);
    xdr.write(&natoms);
    xdr.write(&step);
    xdr.write(&time);
  }


};
