/*
  Simple Grid Class
  (c) 2008 Tod D. Romo,
      Grossfield Lab,
      University of Rochester Medical and Dental School
*/

#if !defined(SGRID_HPP)
#define SGRID_HPP

#include <cmath>
#include <stdexcept>
#include <utility>
#include <boost/iterator/iterator_facade.hpp>

#include <stdexcept>
#include <boost/iterator/iterator_facade.hpp>

#include <loos.hpp>
#include <Coord.hpp>

#include <smetad.hpp>

namespace lab {


  typedef loos::Coord<int> SGridpoint;


  template<class T> class SGrid;

  //! Encapsulates a j-row from an SGrid
  /**
   * This class allows you to access individual columns from the row
   * via operator[].  Only meant to be used in conjunction with
   * SGridPlane.
   */
  template<class T>
  class SGridRow {
  public:
    SGridRow(const long i, SGrid<T>& g) : idx(i), grid(g) { }
  
    T& operator[](const int i) {
      assert(i >= 0 && i < grid.dims[0]);
      assert(grid.ptr != 0);
      return(grid.ptr[idx + i]);
    }

  private:
    long idx;
    SGrid<T>& grid;
  };


  //! Encapsulates an i,j-plane from an SGrid
  /**
   * This class indexes an i,j-plane by row, returning an SGridRow
   * object.  This allows access to elements in the plane through
   * using two operator[] calls, i.e.  a_grid[row][col]
   */
  template<class T>
  class SGridPlane {
  public:
    SGridPlane(const long i, SGrid<T>& g) : idx(i), grid(g) { }

    SGridRow<T> operator[](const int j) {
      assert(j >= 0 && j < grid.dims[1]);
      assert(grid.ptr != 0);
      return(SGridRow<T>(idx + j*grid.dims[0], grid));
    }

  private:
    long idx;
    SGrid<T>& grid;
  };


  //! Random access iterator using the BOOST facade
  /**
   * This iterator generates a family of iterators for an SGrid
   * object.  The iterators are random access and provide a "world()"
   * public method that can return the world coordinates of the
   * current iterator location and a "grid()" that returns the grid
   * coordinates.  Note that this is a method of the iterator.
   * For example,
\code
SGridIterator<double, double> i = grid.begin();
cout << *i << endl;
cout << i.world() << endl;
\endcode
   * the first line will write out the grid value at the origin while
   * the second line will write out the world coords for the origin.
   * On the other hand,
\code
cout << i[0] << endl;
cout << i[0].world() << endl;
\endcode
   * the first line will again write out the grid value at the
   * original, but the second line will be an error...
   */

  // T is the fundamental type (i.e. the type stored in the grid)
  // R is the dereference type (this allows us to use a template to
  // handle both const and non-const iterators).
  //
  // Internally, the iterator position is stored as an index into the
  // block of data managed by the SGrid, and is a long, hence the
  // difference type below...
  template<typename T, typename R>
  class SGridIterator : public boost::iterator_facade<
    SGridIterator<T, R>,
    T,
    boost::random_access_traversal_tag,
    R&,
    long
    >
  {
  public:
    SGridIterator() : src(0), offset(0) { }
    explicit SGridIterator(const SGrid<T>& g, long l) : src(&g), offset(l) { }

    // Converts from any dereference type...
    template<typename S> SGridIterator(const SGridIterator<T, S>& o) : src(o.src), offset(o.offset) { }

    loos::GCoord world() const { return(src->gridToWorld(src->indexToGrid(offset))); }
    loos::GCoord coords() const { return(world()); }
    SGridpoint grid() const { return(src->indexToGrid(offset)); }

  private:

    friend class boost::iterator_core_access;

    // This lets any other SGridIterator access our internal bits...
    template<typename,typename> friend class SGridIterator;

    void increment() { ++offset; }
    void decrement() { --offset; }
    void advance(const long n) { offset += n; }

    template<typename S>
    bool equal(const SGridIterator<T, S>& other) const {
      return(src == other.src && offset == other.offset);
    }

    R& dereference() const {
      if (offset < 0 || offset >= src->dimabc)
        throw(std::range_error("Index out of bounds"));
      return(src->ptr[offset]);
    }

    long distance_to(const SGridIterator<T, R>& other) const {
      if (src != other.src)
        throw(std::logic_error("Grid mismatch"));
      return(other.offset - offset);
    }

    const SGrid<T>* src;
    long offset;
  };


  //! A simple 3D grid class of arbitrary types
  /**
   * This class represents a simple 3D grid of elements.  The grid is
   * located somewhere in realspace and supports conversion from
   * grid-space to real-space and back.  Individual grid elements can
   * be indexed any number of ways.  The SGridPoint type is really a
   * LOOS Coord, but with integer elements.  You can index a
   * grid-point by using an SGridPoint, i.e.
\verbatim
SGridPoint p(1,2,3);
T value = a_grid(p);
\endverbatim
   * SGrid's operator() will automatically make an SGridPoint for you
   * if you pass it multiple integers, i.e.
\verbatim
T value = a_grid(3,2,1);
\endverbatim
   * Note that in this case, the order of the indices is reversed from
   * the SGridPoint constructor.  This may seem a little confusing at
   * first, but it was meant to be consistent with how you would
   * access the grid with the operator[] (shown below).
   * As an alternative, you can access the grid as though it was a linear
   * array by passing a long to operator(), i.e.
\verbatim
T value = a_grid(1024);
\endverbatim
   * Finally, you can slice up a grid using operator[], so the first
   * example above could be rewritten as:
\verbatim
T value = a_grid[3][2][1];
\endverbatim
   * This method is not terribly efficient, unless you are going to
   * operate say just a plane.
   *
   * Storing and loading grids is very easy.  Just use the << and >>
   * operators...
   */
  template<class T>
  class SGrid {
  public:

    typedef T                          value_type;
    typedef SGridIterator<T, T>        iterator;
    typedef SGridIterator<T, const T>  const_iterator;

    friend class SGridRow<T>;
    friend class SGridPlane<T>;

    friend class SGridIterator<T, T>;
    friend class SGridIterator<T, const T>;

    typedef std::pair<int, int>        Range;

    //! Empty grid
    SGrid() : ptr(0), _gridmin(loos::GCoord(0,0,0)), _gridmax(loos::GCoord(0,0,0)),
	      dims(SGridpoint(0,0,0)) { }

    //! Create a grid with explicit location in realspace and dimensions
    /**
     * The passed GCoords are the corners of the grid in real-space.
     * The SGridpoint defines the size of each dimension of the grid.
     */
    SGrid(const loos::GCoord& gmin, const loos::GCoord& gmax, const SGridpoint& griddims) :
      ptr(0), _gridmin(gmin), _gridmax(gmax), dims(griddims) { init(); }

    //! Creates a grid with an explicit location and uniform size in
    //! all dimensions
    SGrid(const loos::GCoord& gmin, const loos::GCoord& gmax, const int dim) :
      ptr(0), _gridmin(gmin), _gridmax(gmax), dims(dim, dim, dim) { init(); }

    //! Copies an existing grid
    SGrid(const SGrid<T>& g) : ptr(0), _gridmin(g._gridmin), _gridmax(g._gridmax),
			       dims(g.dims) {
      init();
      memcpy(ptr, g.ptr, dimabc * sizeof(T));
      meta_ = g.meta_;
    }

    void resize(const loos::GCoord& gmin, const loos::GCoord& gmax, const SGridpoint& griddims) {
      if (ptr != 0)
        delete[] ptr;
      _gridmin = gmin;
      _gridmax = gmax;
      dims = griddims;
      init();
    }

    //! This is a "deep" copy of grid
    const SGrid<T>& operator=(const SGrid<T>& g) {
      if (this == &g)
	return(*this);

      if (ptr) {
	delete[] ptr;
	ptr = 0;
      }
      _gridmin = g._gridmin;
      _gridmax = g._gridmax;
      dims = g.dims;
      init();
      if (dimabc != 0)
        memcpy(ptr, g.ptr, dimabc * sizeof(T));

      meta_ = g.meta_;

      return(*this);
    }


    ~SGrid() { delete[] ptr; }



    SGrid<T> subset(const Range& c, const Range& b, const Range& a) {
      SGridpoint dim;

      dim.x(a.second - a.first + 1);
      dim.y(b.second - b.first + 1);
      dim.z(c.second - c.first + 1);

      loos::GCoord bottom = gridToWorld(SGridpoint(a.first, b.first, c.first));
      loos::GCoord top = gridToWorld(SGridpoint(a.second, b.second, c.second));

      SGrid<T> sub(bottom, top, dim);
      for (int k=0; k<dim.z(); ++k)
        for (int j=0; j<dim.y(); ++j)
          for (int i=0; i<dim.x(); ++i)
            sub(k, j, i) = operator()(k+c.first, j+b.first, i+a.first);

      return(sub);
    }

    //! Zero out all elements
    void zero(void) { for (long i=0; i<dimabc; i++) ptr[i] = 0; }
  

    //! Takes an SGridPoint and returns the "linear" index into the
    //! grid
    long gridToIndex(const SGridpoint v) const {
      return( (v[2] * dims[1] + v[1]) * dims[0] + v[0] );
    }


    //! Converts a real-space coordinate into grid coords
    SGridpoint gridpoint(const loos::GCoord& x) const {
      SGridpoint v;

      for (int i=0; i<3; i++) {
	long k = static_cast<long>(floor( (x[i] - _gridmin[i]) * delta[i] + 0.5 ));
	v[i] = k;
      }

      return(v);
    }
    
    //! Converts a real-space coordinate into grid coords
    SGridpoint gridpoint(const double z, const double y, const double x) const {
      loos::GCoord c(x,y,z);
      return(gridpoint(c));
    }


    //! Checks to make sure the gridpoint lies within the grid boundaries
    bool inRange(const SGridpoint& g) const {
      for (int i=0; i<3; i++)
	if (g[i] < 0 || g[i] >= dims[i])
	  return(false);

      return(true);
    }

    bool inRange(const int k, const int j, const int i) const {
      return(inRange(SGridpoint(i, j, k)));
    }

    //! Access the grid element indexed by k, j, i

    // Bypassing the Coord<>::operator[] saves us a considerable
    // amount of time...hence the Coord<>::x(), etc...
     
    T& operator()(const int k, const int j, const int i) {
      long x = ((k * dims.y()) + j) * dims.x() + i;
      assert(x < dimabc);
      return(ptr[x]);
    }

    //! Access the grid element indexed by the SGridPoint
    T& operator()(const SGridpoint& v) {
      for (int i=0; i<3; i++)
	assert(v[i] >= 0 && v[i] < dims[i]);
      return(ptr[ (v[2] * dims.y() + v[1]) * dims.x() + v[0] ]);
    }

    //! Access the element indexed by i, assuming the grid to be a big
    //! linear array
    T& operator()(const long i) {
      assert(i >= 0 && i < dimabc);
      return(ptr[i]);
    }

    //! Converts \a x into grid coords, then accesses that element
    T& operator()(const loos::GCoord& x) {
      return(operator()(gridpoint(x)));
    }


    // Const versions...

      
    const T& operator()(const int k, const int j, const int i) const {
      long x = ((k * dims.y()) + j) * dims.x() + i;
      assert(x < dimabc);
      return(ptr[x]);
    }

    const T& operator()(const SGridpoint& v) const {
      for (int i=0; i<3; i++)
	assert(v[i] >= 0 && v[i] < dims[i]);
      return(ptr[ (v[2] * dims.y() + v[1]) * dims.x() + v[0] ]);
    }

    const T& operator()(const long i) const {
      assert(i >= 0 && i < dimabc);
      return(ptr[i]);
    }


    const T& operator()(const loos::GCoord& x) const {
      return(operator()(gridpoint(x)));
    }


    // Slicin' und dicin'...

    //! Returns the kth plane from the grid
    SGridPlane<T> operator[](const int k) {
      assert(k >= 0 && k < dims[2]);
      return(SGridPlane<T>(k * dimab, *this));
    }


    //! Converts grid coords to real-space (world) coords
    loos::GCoord gridToWorld(const SGridpoint& v) const {
      loos::GCoord c;
    
      for (int i=0; i<3; i++)
	c[i] = static_cast<loos::greal>(v[i]) / delta[i] + _gridmin[i];

      return(c);
    }

    //! Calculates the grid coords from a linear index
    SGridpoint indexToGrid(const long idx) const {
      int c = idx / dimab;
      int r = idx % dimab;
      int b = r / dims[0];
      int a = r % dims[0];

      return(SGridpoint(a, b, c));
    }
  

    //! Squared-distance (in real-space) between two grid coords
    double gridDist2(const SGridpoint& u, const SGridpoint& v) {
      loos::GCoord x = gridToWorld(u);
      loos::GCoord y = gridToWorld(v);

      return(x.distance2(y));
    }

    //! Linear distance (in real-space) between two grid coords
    double gridDist(const SGridpoint& u, const SGridpoint& v) {

      return(sqrt(gridDist2(u, v)));
    }

    //! Returns a vector of SGridPoints that lie within a box
    //! containing a sphere of radius \a r
    std::vector<SGridpoint> withinBoxRadius(const double r, const loos::GCoord& u, const int pad = 0) const {

      SGridpoint a = gridpoint(u - r);
      SGridpoint b = gridpoint(u + r);

      std::vector<SGridpoint> res;
      for (int k=a[2]-pad; k <= b[2]+pad; k++) {
	if (k < 0 || k >= dims[2])
	  continue;
	for (int j=a[1]-pad; j <= b[1]+pad; j++) {
	  if (j < 0 || j >= dims[1])
	    continue;
	  for (int i=a[0]-pad;i <= b[0]+pad; i++) {
	    if (i < 0 || i >= dims[0])
	      continue;
	    res.push_back(SGridpoint(i, j, k));
	  }
	}
      }
    
      return(res);
    }


    //! Returns a vector of SGridPoints that lie within a given
    //! radius of a real-space coord
    /**
     * This first converts the real-space coords into grid-space, then
     * finds the bounding box for the sphere with radius \r.  It then
     * checks each grid-point to see if it actually lies on the
     * surface or within the sphere and returns a vector of
     * SGridPoints that do.
     */
    std::vector<SGridpoint> withinRadius(const double r, const loos::GCoord& u) const {
      std::vector<SGridpoint> pts = withinBoxRadius(r, u);
      std::vector<SGridpoint> res;
      std::vector<SGridpoint>::iterator i;

      double r2 = r*r;

      for (i=pts.begin(); i != pts.end(); i++) {
	loos::GCoord v = gridToWorld(*i);
	if (u.distance2(v) <= r2)
	  res.push_back(*i);
      }

      return(res);
    }


    //! Applies a functor to all grid points that lie on or within a
    //! sphere of radius \a r about the real-space coords \a u.
    /*
     * Here, a functor \a f is called for all grid points that lie on
     * or within a sphere of radius \a r about the real-space coords
     * \a u.  The functor take two parameters: the value at the gid
     * point and the distance from that grid point to the center of
     * the sphere (in real-space units)
     */
    template<typename Func>
    void applyWithinRadius(const double r, const loos::GCoord& u, const Func& f) {
      SGridpoint a = gridpoint(u - r);
      SGridpoint b = gridpoint(u + r);
      double r2 = r * r;

      std::vector<SGridpoint> res;
      for (int k=a[2]; k <= b[2]; k++) {
	if (k < 0 || k >= dims[2])
	  continue;
	for (int j=a[1]; j <= b[1]; j++) {
	  if (j < 0 || j >= dims[1])
	    continue;
	  for (int i=a[0];i <= b[0]; i++) {
	    if (i < 0 || i >= dims[0])
	      continue;
	    
	    SGridpoint point(i, j, k);
	    loos::GCoord v = gridToWorld(point);
	    double d = u.distance2(v);
	    if (d <= r2)
	      f(operator()(v), d);
	  }
	}
      }


    }


    void scale(const T val) {
      for (long i = 0; i < dimabc; ++i)
        ptr[i] *= val;
    }

    void clear(const T val = 0) {
      for (long i = 0; i < dimabc; ++i)
        ptr[i] = val;
    }


    SGridpoint gridDims(void) const { return(dims); }
    loos::GCoord minCoord(void) const { return(_gridmin); }
    loos::GCoord maxCoord(void) const { return(_gridmax); }
    loos::GCoord gridDelta(void) const { return(delta); }

    long maxGridIndex(void) const { return(dimabc); }
    long size() const { return(dimabc); }

    //! Write out a grid
    /**
     * This simply linearizes the grid and writes it out (with
     * associated metadata).  This means that if your grid is mostly
     * empty, it will be rather wasteful.  On the other handle, it IS
     * a simple grid implementation...
     */
    friend std::ostream& operator<<(std::ostream& os, const SGrid<T>& grid) {
      os << "# SGrid-1.1\n";
      os << grid.meta_;
      os << grid.dims << std::endl;
      os << grid._gridmin << std::endl;
      os << grid._gridmax << std::endl;
      //      os << boost::format("(%.8f,%.8f,%.8f)\n") % grid._gridmin[0] % grid._gridmin[1] % grid._gridmin[2];
      //      os << boost::format("(%.8f,%.8f,%.8f)\n") % grid._gridmax[0] % grid._gridmax[1] % grid._gridmax[2];

      return(os.write(reinterpret_cast<char*>(grid.ptr), sizeof(T) * grid.dimabc));
    }

    //! Read in a grid
    /**
     * Any existing grid will get clobbered--replaced by the grid
     * being read in.
     */
    friend std::istream& operator>>(std::istream& is, SGrid<T>& grid) {
      std::string s;

      std::getline(is, s);
      if (s != "# SGrid-1.1")
	throw(std::runtime_error("Bad input format for SGrid  - " + s));

      is >> grid.meta_;

      if (grid.ptr)
	delete[] grid.ptr;

      is >> grid.dims;
      is >> grid._gridmin;
      is >> grid._gridmax;
      if (is.get() != '\n')  // Pull trailing newline off of input...
        throw(std::runtime_error("Grid parse error in header"));

      grid.init();
      is.read(reinterpret_cast<char*>(grid.ptr), sizeof(T) * grid.dimabc);
      if (is.fail() || is.eof())
        throw(std::runtime_error("Grid read error"));

      return(is);
    }

    iterator begin() { return(iterator(*this, 0)); }
    iterator end() { return(iterator(*this, dimabc)); }

    const_iterator begin() const { return(const_iterator(*this, 0)); }
    const_iterator end() const { return(const_iterator(*this, dimabc)); }

    bool empty() const { return(dimabc == 0); }

    void setMetadata(const std::string& s) { meta_.set(s); }
    void addMetadata(const std::string& s) { meta_.add(s); }

    SMetaData metadata() const { return(meta_); }
    void metadata(const SMetaData& m) { meta_ = m; }
    

  private:
    void init(void) {
      dimab = dims[0]*dims[1];
      dimabc = dimab * dims[2];

      for (int i=0; i<3; i++)
	delta[i] = (dims[i] - 1)/ (_gridmax[i] - _gridmin[i]);

      if (dimabc != 0) {
        ptr = new T[dimabc];
        zero();
      } else
        ptr = 0;
    }

  private:
    T* ptr;
    loos::GCoord _gridmin, _gridmax, delta;
    SGridpoint dims;
    long dimabc, dimab;

    SMetaData meta_;
  };


};

#endif
