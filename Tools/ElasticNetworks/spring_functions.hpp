/*
  Spring Functions from Elastic Network Tools in LOOS
  (c) 2010 Tod D. Romo, Grossfield Lab, URMC
*/



#if !defined(SPRING_FUNCTIONS_HPP)
#define SPRING_FUNCTIONS_HPP




#include <loos.hpp>



// These classes define the various possible spring functions used in
// creating the Hessian.  All derived from the SpringFunction base
// class.  This class returns a 3x3 DoubleMatrix containing the spring
// constants...
//
// The SpringFunction::constant() function takes the coords of the two
// nodes plus their difference vector (since it'll almost always be
// computed prior to calling SpringFunction, no since in recomputing
// it).
//
//
// * Misc notes *
//
// The setConstants() function is destructive to its parameter.  I'm
// not sure if this is the best way to implement this, so don't depend
// on it staying this way.  The vector is treated as a LIFO stack.  It
// is expected that where there are multiple constants, they will be
// pushed in the order they would have been passed had setConstants()
// taken a regular argument list, i.e.
//    head --> rcut, k1, k2, k3, k4  <-- tail


class SpringFunction {
public:
  SpringFunction() : warned(false) { }
  virtual ~SpringFunction() { }
  virtual std::string name() const =0;

  // Sets the internal spring constants (destructively).  Returns true
  // if constants are "valid"
  virtual bool setConstants(std::vector<double>& konst) =0;
  // Returns the # of spring constants used by this spring function
  virtual uint numberOfConstants() const =0;

  
  
  virtual loos::DoubleMatrix constant(const loos::GCoord& u, const loos::GCoord& v, const loos::GCoord& d) =0;

protected:

  // Check for negative spring-constants.  If found, issue a one-time
  // warning, but only for this specific spring function...

  double checkConstant(double d) {
    if (d < 0.0 && !warned) {
      warned = true;
      std::cerr << "Warning- negative spring constants found in " << name() << ".  Setting to 0.\n";
      d = 0.0;
    }

    return(d);
  }

private:
  bool warned;
};




// Many spring functions are actually uniform over the 3x3 matrix.
// The UniformSpringFunction uses the NVI idiom
// (http://www.gotw.ca/publications/mill18.htm) to allow subclasses to
// only return a double value, which then gets copied into all
// elements of the 3x3 matrix.
//
// Note: this means you override the constantImpl() implementation
// function, NOT the public constant() function.

class UniformSpringFunction : public SpringFunction {
public:

  loos::DoubleMatrix constant(const loos::GCoord& u, const loos::GCoord& v, const loos::GCoord& d) {
    double k = checkConstant(constantImpl(u, v, d));
    loos::DoubleMatrix B(3, 3);
    for (uint i=0; i<9; ++i)
      B[i] = k;

    return(B);
  }

private:
  virtual double constantImpl(const loos::GCoord& u, const loos::GCoord& v, const loos::GCoord& d) =0;
};



// The following are all uniform spring constants...


// Basic distance cutoff for "traditional" ENM

class DistanceCutoff : public UniformSpringFunction {
public:
  DistanceCutoff(const double& r) : radius(r*r) { }
  DistanceCutoff() : radius(15.0*15.0) { }

  std::string name() const { return("DistanceCutoff"); }

  bool setConstants(std::vector<double>& konst) {
    radius = konst.back();
    radius *= radius;
    konst.pop_back();
    return(radius != 0.0);
  }

  uint numberOfConstants() const { return(1); }


  double constantImpl(const loos::GCoord& u, const loos::GCoord& v, const loos::GCoord& d) {
    double s = d.length2();
    if (s <= radius)
      return(1./s);

    return(0.0);
  }

private:
  double radius;
};


// Distance weighting (i.e. ||u-v||^p)

class DistanceWeight : public UniformSpringFunction {
public:
  DistanceWeight(const double p) : power(p) { }
  DistanceWeight() : power(-2.0) { }

  std::string name() const { return("DistanceWeight"); }


  bool setConstants(std::vector<double>& konst) {
    power = konst.back();
    konst.pop_back();
    return(power < 0.0);
  }

  uint numberOfConstants() const { return(1); }


  double constantImpl(const loos::GCoord& u, const loos::GCoord& v, const loos::GCoord& d) {
    double s = d.length();
    return(pow(s, power));
  }

private:
  double power;
};




// Exponential distance weighting (i.e. exp(k * ||u-v||))

class ExponentialDistance : public UniformSpringFunction {
public:
  ExponentialDistance(const double s) : scale(s) { }
  ExponentialDistance() : scale(-2.0) { }

  std::string name() const { return("ExponentialDistance"); }

  bool setConstants(std::vector<double>& konst) {
    scale = konst.back();
    konst.pop_back();
    return(scale != 0.0);
  }

  uint numberOfConstants() const { return(1); }


  double constantImpl(const loos::GCoord& u, const loos::GCoord& v, const loos::GCoord& d) {
    double s = d.length();
    return(exp(scale * s));
  }

private:
  double scale;
};





// HCA method (see Hinsen et al, Chem Phys (2000) 261:25-37)
// Note: The defaults are the original Hinsen constants...

class HCA : public UniformSpringFunction {
public:
  HCA(const double rc, const double a, const double b, const double c, const double d) :
    rcut(rc), k1(a), k2(b), k3(c), k4(d) { }

  HCA() :
    rcut(4.0), k1(205.5), k2(571.2), k3(305.9e3), k4(6.0) { }

  std::string name() const { return("HCA"); }

  bool setConstants(std::vector<double>& konst) {
    k4 = konst.back();
    konst.pop_back();

    k3 = konst.back();
    konst.pop_back();

    k2 = konst.back();
    konst.pop_back();

    k1 = konst.back();
    konst.pop_back();

    rcut = konst.back();
    konst.pop_back();

    return( rcut >= 0.0 && k4 >= 0.0 );
  }

  uint numberOfConstants() const { return(5); }


  double constantImpl(const loos::GCoord& u, const loos::GCoord& v, const loos::GCoord& d) {
    double s = d.length();
    
    if (s <= rcut)
      return(k1 * s - k2);

    return(k3 * pow(s, -k4));
  }

private:
  double rcut, k1, k2, k3, k4;
};




SpringFunction* springFactory(const std::string& spring_desc);


#endif
