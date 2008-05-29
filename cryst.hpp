/*
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Crystallographic unit cell params...
*/



#if !defined(CRYST_HPP)
#define CRYST_HPP

class UnitCell {
public:
  UnitCell() : _a(1.0), _b(1.0), _c(1.0), _alpha(90.0), _beta(90.0), _gamma(90.0), sgroup("P1"), zval(1) { }

  greal a(void) const { return(_a); }
  void a(const greal x) { _a = x; }

  greal b(void) const { return(_b); }
  void b(const greal x) { _b = x; }

  greal c(void) const { return(_c); }
  void c(const greal x) { _c = x; }

  greal alpha(void) const { return(_alpha); }
  void alpha(const greal x) { _alpha = x; }

  greal beta(void) const { return(_beta); }
  void beta(const greal x) { _beta = x; }

  greal gamma(void) const { return(_gamma); }
  void gamma(const greal x) { _gamma = x; }

  string spaceGroup(void) const { return(sgroup); }
  void spaceGroup(const string s) { sgroup = s; }

  int z(void) const { return(zval); }
  void z(const int i) { zval = i; }

  friend ostream& operator<<(ostream& os, const UnitCell& u) {
    os << "<UNITCELL A='" << u._a << "' B='" << u._b << "' C='" << u._c << "' ALPHA='";
    os << u._alpha << "' BETA='" << u._beta << "' GAMMA='" << u._gamma << "' SPACEGROUP='";
    os << u.sgroup << "' Z='" << u.zval << "'/>";
    return(os);
  }



private:
  greal _a, _b, _c, _alpha, _beta, _gamma;
  string sgroup;
  int zval;
};


#endif
