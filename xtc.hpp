#if !defined(XTC_HPP)
#define XTC_HPP
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

  class XTC : public Trajectory {

    struct Header {
      uint magic;
      uint natoms, step;
      float time, box[9];
    };

    static const int magicints[];
    static const int firstidx, lastidx;
    typedef float    xtc_t;

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

    void updateGroupCoords(AtomicGroup& g);

    std::vector<GCoord> coords(void) { return(coords_); }

  private:

    internal::XDR xdr_file;
    std::vector<size_t> frame_indices;
    uint natoms_;
    GCoord box;
    double precision_;
    std::vector<GCoord> coords_;


  private:

    bool parseFrame(void);
    int sizeofint(int);
    int sizeofints(uint*, const uint);
    int decodebits(int*, uint);
    void decodeints(int*, const int, int, uint*, int*);
    bool readFrameHeader(Header&);
    void scanFrames(void);
    
    void seekNextFrameImpl(void) { }
    void seekFrameImpl(uint);
    void rewindImpl(void) { ifs()->clear(); ifs()->seekg(0); }
    bool readCompressedCoords(std::vector<GCoord>&);

  };

}



#endif
