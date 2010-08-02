#if !defined(SPRING_FUNCTIONS_HPP)
#define SPRING_FUNCTIONS_HPP




#include <loos.hpp>



class SpringFunction {
public:
  SpringFunction() : warned(false) { }
  virtual ~SpringFunction() { }
  virtual std::string name() const =0;
  
  
  virtual loos::DoubleMatrix constant(const loos::GCoord& u, const loos::GCoord& v, const loos::GCoord& d) =0;

protected:
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




class DistanceCutoff : public UniformSpringFunction {
public:
  DistanceCutoff(const double& r) : radius(r*r) { }

  std::string name() const { return("DistanceCutoff"); }

  double constantImpl(const loos::GCoord& u, const loos::GCoord& v, const loos::GCoord& d) {
    double s = d.length2();
    if (s <= radius)
      return(1./s);

    return(0.0);
  }

private:
  double radius;
};



class DistanceWeight : public UniformSpringFunction {
public:
  DistanceWeight(const double p) : power(p) { }

  std::string name() const { return("DistanceWeight"); }

  double constantImpl(const loos::GCoord& u, const loos::GCoord& v, const loos::GCoord& d) {
    double s = d.length();
    return(pow(s, power));
  }

private:
  double power;
};



class ExponentialDistance : public UniformSpringFunction {
public:
  ExponentialDistance(const double s) : scale(s) { }

  std::string name() const { return("ExponentialDistance"); }

  double constantImpl(const loos::GCoord& u, const loos::GCoord& v, const loos::GCoord& d) {
    double s = d.length();
    return(exp(scale * s));
  }

private:
  double scale;
};



class HCA : public UniformSpringFunction {
public:
  HCA(const double rc, const double a, const double b, const double c, const double d) :
    rcut(rc), k1(a), k2(b), k3(c), k4(d) { }

  HCA() :
    rcut(4.0), k1(205.5), k2(571.2), k3(305.9e3), k4(6.0) { }

  std::string name() const { return("HCA"); }

  double constantImpl(const loos::GCoord& u, const loos::GCoord& v, const loos::GCoord& d) {
    double s = d.length();
    
    if (s <= rcut)
      return(k1 * s - k2);

    return(k3 * pow(s, -k4));
  }

private:
  double rcut, k1, k2, k3, k4;
};



#endif
