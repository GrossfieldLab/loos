#include <iostream>
#include <stdexcept>



#include "Coord.hpp"
#include "Geometry.hpp"

#define TYPE float



int main(int argc, char *argv[]) {


  Coord<TYPE> u(3, 4, 5), v(1, -1, 0);

  cout << "U: " << u << ", V: " << v << endl;
  cout << "u+10 : 10+u\n";

  cout << u+10 << endl;
  cout << 10+u << endl;

  cout << "-- const tests\n";
  cout << "u-10 : 10-u : -u\n";
  
  cout << u-10 << endl;
  cout << 10-u << endl;
  cout << -u << endl;


  cout << "-- add/sub tests\n";

  cout << "u-v : v-u : u+v\n";
  cout << u-v << endl;
  cout << v-u << endl;
  cout << u+v << endl;

  cout << "-- mult tests\n";
  cout << "u*2 : u/2\n";

  cout << u*2 << endl;
  cout << u/2 << endl;

  cout << "-- dot and cross tests\n";
  Coord<TYPE> a(1,2,3), b(4, 5, 6);
  cout << "A: " << a << ", B: " << b << endl;
  cout << "(dot product) : (cross product)\n";
  cout << a * b << endl;
  cout << (a ^ b) << endl;

  cout << "-- distance tests\n";

  cout << "d^2(a->b) = " << a.distance2(b) << endl;

  cout << "-- Mod tests";
  Coord<TYPE> x(2,2,4);
  cout << "X: " << x << endl;
  cout << "b % 2 : b % x\n";
  cout << b % 2 << endl;
  cout << b % x << endl;

  cout << "------\n";
  cout << "normalized_b = " << b / b.length() << endl;

  cout << "-- Box tests\n";
  u.set(2, 0, 0);
  v.set(8, 0, 0);
  Coord<TYPE> box(10,10,10);

  cout << u << " -> " << v << " in box " << box << " = " << u.distance(v, box) << "\n";

  
  Coord<double> a1(1,0,0), a2(0,0,0), a3(0,1,0), a4(0,1,1);
  cout << angle(a1,a2,a3) << endl;
  cout << angle(a1,a2,a4) << endl;
  cout << angle(a1,a3,a4) << endl;
  cout << torsion(a1,a2,a3,a4) << endl;


}
