/*
  Xplor EDM writer (writing a grid)
  (c) 2008 Tod D. Romo,
      Grossfield Lab,
      University of Rochester Medical and Dental School
*/


#if !defined(XPLOREDMWRITER_HPP)
#define XPLOREDMWRITER_HPP

#include <loos.hpp>

#include <DensityGrid.hpp>

namespace loos {

  //! Functor for writing out ASCII formatted X-Plor electron density maps
  template<class T>
  struct XEDMWriter {
    XEDMWriter(std::ostream& os) : i(0), mos(os), fmt(5)
    {
      fmt.scientific().width(12).right();
    }

    void operator()(const T d) {
      mos << fmt(d);
      if (i++ >= 5) {
	mos << std::endl;
	i = 0;
      }
    }

    //! Special m.f. for starting a new frame
    void frame(const int k) {
      if (i != 0) {
	i = 0;
	mos << std::endl;
      }
      mos << std::setw(8) << k << std::endl;
    }


    int i;
    std::ostream& mos;
    loos::Fmt fmt;
  };



  //! Write out an SGrid as an ASCII formatted X-PLOR electron density map
  template<class T> void writeXplorEDM(std::ostream& os, SGrid<T>& grid) {
    loos::GCoord gridmin = grid.minCoord();
    loos::GCoord gridmax = grid.maxCoord();
    loos::GCoord delta = grid.gridDelta();
    loos::GCoord gridsize;
    SGridpoint dims = grid.gridDims();

    SGridpoint mins, maxs, nas;

    // Handle the header...
    for (int i=0; i<3; i++) {
      mins[i] = static_cast<int>(floor(gridmin[i] * delta[i]));
      maxs[i] = static_cast<int>(floor(gridmax[i] * delta[i]));
      gridsize[i] = dims[i] / delta[i];
    }
    nas = dims;

//     os << std::endl << std::setw(8) << 1 << std::endl;
//     os << "XPLOR-EDM FROM SGRID\n";

    // Special handling for grid meta-data...
    SMetaData meta = grid.metadata();
    os << std::endl << std::setw(8) << meta.size() << std::endl;
    for (SMetaData::iterator i = meta.begin(); i != meta.end(); ++i)
      os << *i << std::endl;

    for (int i=0; i<3; i++)
      os << std::setw(8) << nas[i] << std::setw(8) << mins[i] << std::setw(8) << maxs[i];
    os << std::endl;

    loos::Fmt fc(5);
    fc.width(12).scientific();

    // Assume our "crystal" is orthonormal...
    os << fc(gridsize[0]) << fc(gridsize[1]) << fc(gridsize[2]) << fc(90.0) << fc(90.0) << fc(90.0) << std::endl;
    os << "ZYX\n";

    // Instantiate the writing functor...
    XEDMWriter<T> writer(os);

    // The format writes out the map a plane at a time, so we extract
    // a plane via operator[] and operate on that...

    for (int k=0; k<dims[2]; k++) {
      SGridPlane<T> plane = grid[k];

      // Prime the output
      writer.frame(k);

      for (int j=0; j<dims[1]; j++)
	for (int i=0; i<dims[0]; i++)
	  writer(plane[j][i]);
    }

    os << std::endl << std::endl;
  }

};


#endif
