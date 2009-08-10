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


#if !defined(TRR_HPP)
#define TRR_HPP
#include <ios>
#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

#include <loos_defs.hpp>

#include <Coord.hpp>
#include <xdr.hpp>
#include <AtomicGroup.hpp>
#include <Trajectory.hpp>

#include <boost/format.hpp>


namespace loos {

  class TRR : public Trajectory {
    static const int magic;
    static const int DIM;

    struct Header {
      bool bDouble;
      int ir_size;
      int e_size;
      int box_size;
      int vir_size;
      int pres_size;
      int top_size;
      int sym_size;
      int x_size;
      int v_size;
      int f_size;

      int natoms;
      int step;
      int nre;
      float tf;
      float lambdaf;
      double td;
      double lambdad;
    };

  public:
    explicit TRR(const std::string& s) : Trajectory(s), xdr_file(ifs()) {
      init();
      rewindImpl();
      parseFrame();
    }

    explicit TRR(const char* p) : Trajectory(p), xdr_file(ifs()) {
      init();
      rewindImpl();
      parseFrame();
    }
    explicit TRR(std::iostream& is) : Trajectory(is), xdr_file(ifs()) {
      init();
      rewindImpl();
      parseFrame();
    }

    uint natoms(void) const { return(hdr_.natoms); }
    float timestep(void) const { return(0); }
    uint nframes(void) const { return(frame_indices.size()); }
    bool hasPeriodicBox(void) const { return(hdr_.box_size != 0); }
    GCoord periodicBox(void) const { return(box); }
    
    void updateGroupCoords(AtomicGroup& g);
    
    std::vector<GCoord> coords(void) { return(coords_); }

    // TRR specific attributes...
    std::vector<double> virial(void) const { return(vir_); }
    std::vector<double> pressure(void) const { return(pres_); }
    std::vector<GCoord> velocities(void) const { return(velo_); }
    std::vector<GCoord> forces(void) const { return(forc_); }

    bool isDouble(void) const { return(hdr_.bDouble); }
    bool hasVirial(void) const { return(hdr_.vir_size != 0); }
    bool hasPressure(void) const { return(hdr_.pres_size != 0); }
    bool hasCoords(void) const { return(hdr_.x_size != 0); }
    bool hasVelocities(void) const { return(hdr_.v_size != 0); }
    bool hasForces(void) const { return(hdr_.f_size != 0); }

    double time(void) const { return( hdr_.bDouble ? hdr_.td : hdr_.tf ); }
    double lambda(void) const { return( hdr_.bDouble ? hdr_.lambdad : hdr_.lambdaf); }

    int step(void) const { return(hdr_.step); }


  private:
    void init(void);
    int floatSize(Header& h);
    bool readHeader(Header& h);

    template<typename T>
    void readBlock(std::vector<double>& v, const uint n, const std::string& msg) {
      T buf[n];
      uint i = xdr_file.read(buf, n);
      if (i != n)
        throw(std::runtime_error("Unable to read TRR " + msg));
      for (i=0; i<n; ++i)
        v.push_back(buf[i]);
    }

    // By default, this converts to angstroms...
    template<typename T>
    void readBlock(std::vector<GCoord>& v, const uint n, const std::string& msg) {
      T* buf = new T[n];
      if (buf == 0)
        throw(std::runtime_error("Out of memory"));

      if (xdr_file.read(buf, n) != n) {
        delete[] buf;
        throw(std::runtime_error("Unable to read TRR " + msg));
      }
      for (uint i=0; i<n; i += DIM)
        v.push_back(GCoord(buf[i], buf[i+1], buf[i+2]) * 10.0);

      delete[] buf;
    }


    template<typename T>
    bool readRawFrame() {
      
      // Should it be an error/warning if the frame has different
      // number of atoms from what was previously read?
      natoms_ = hdr_.natoms;

      // Clear data first...
      box_.clear();
      vir_.clear();
      pres_.clear();
      velo_.clear();
      forc_.clear();
      coords_.clear();

      if (hdr_.box_size) {
        readBlock<T>(box_, DIM*DIM, "box");
        box = GCoord(box_[0], box_[4], box_[8]) * 10.0;   // Convert
                                                          // to angstroms
      }

      if (hdr_.vir_size)
        readBlock<T>(vir_, DIM*DIM, "virial");
      if (hdr_.pres_size)
        readBlock<T>(pres_, DIM*DIM, "pressure");
      
      if (hdr_.x_size)
        readBlock<T>(coords_, hdr_.natoms * DIM, "Coordinates");

      if (hdr_.v_size)
        readBlock<T>(velo_, hdr_.natoms * DIM, "Velocities");
      
      if (hdr_.f_size)
        readBlock<T>(forc_, hdr_.natoms * DIM, "Forces");

      
      return(! ((xdr_file.get())->fail() || (xdr_file.get())->eof()) );
    }

    void rewindImpl(void) { ifs()->clear(); ifs()->seekg(0, std::ios_base::beg); }
    void seekNextFrameImpl(void) { }
    void seekFrameImpl(uint);


    bool parseFrame(void) {
      if (!readHeader(hdr_))
        return(false);

      if (hdr_.bDouble)
        return(readRawFrame<double>());

      return(readRawFrame<float>());
    }


  private:
    internal::XDR xdr_file;
    std::vector<GCoord> coords_;
    GCoord box;
    uint natoms_;
    std::vector<size_t> frame_indices;

    std::vector<double> box_;
    std::vector<double> vir_;
    std::vector<double> pres_;
    std::vector<GCoord> velo_;
    std::vector<GCoord> forc_;

    Header hdr_;
  };



};


#endif

