/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo
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


#include <loos.hpp>


// # of iterations
static const int maxcount = 10000;

// # of pseudo-atoms
static const int nparticles = 100;

// Threshold for coarse tests...  This is comparison between the
// principal axes and the transformed coordinate axes
static const double threshold = 0.2;

// Fine tests threshold...  This is where we calculate PA's for the
// group, then randomly rotate and translate it again and compare the
// new PA's against the transformed old PA's (i.e. the Alan Test)
static const double fine_threshold = 1e-10;

// Size of the axes of the ellipsoid
static const greal ellipsoid[3] = {10.0, 5.0, 20.0};

// Exit the program with an error code if there is a failure in a test.
static const bool exit_on_failure = false;

// Show the raw results..
static const bool show_results = false;


static loos::base_generator_type& rng = loos::rng_singleton();
static const double pi = 4.0*atan2(1,1);


XForm randomXForm(void) {
  XForm M;

  boost::uniform_real<> map(0.0, 1.0);
  boost::variate_generator<loos::base_generator_type&, boost::uniform_real<> > randy(rng, map);

  greal xt = 50 * randy() - 25;
  greal yt = 50 * randy() - 25;
  greal zt = 50 * randy() - 25;

  M.translate(xt, yt, zt);
  M.push();
  M.identity();
  M.rotate('z', 360*randy());
  M.rotate('y', 360*randy());
  M.rotate('z', 360*randy());
  
  return(M);
}


AtomicGroup createGroup(const int natoms, const float a, const float b, const float c) {
  AtomicGroup A;
  XForm M;
  
  boost::uniform_real<> map(0.0, 1.0);
  boost::variate_generator<loos::base_generator_type&, boost::uniform_real<> > randy(rng, map);
  
  int i;
  for (i=0; i<natoms; i++) {
    greal theta = 2*pi*randy();
    greal phi = pi*randy();
    
    greal x = a*cos(theta)*sin(phi);
    greal y = b*sin(theta)*sin(phi);
    greal z = c*cos(phi);
    
      pAtom pa(new Atom(i, "CA", GCoord(x, y, z)));
    A += pa;
  }
  
  return(A);
}

GCoord matchSigns(const GCoord& a, const GCoord& b) {
  GCoord c;
  int i;

  for (i=0; i<3; i++)
    c[i] = (a[i]<0.0 && b[i]>0.0 || a[i]>0.0 && b[i]<0.0) ? -b[i] : b[i];
      
  return(c);
}


vector<GCoord> computeRotation(XForm& M) {
  vector<GCoord> result;

  GCoord v(0, 0, 1);
  v = M.transform(v);
  v /= v.length();
  result.push_back(v);

  v.set(1, 0, 0);
  v = M.transform(v);
  v /= v.length();
  result.push_back(v);

  v.set(0, 1, 0);
  v = M.transform(v);
  v /= v.length();
  result.push_back(v);

  return(result);
}



vector<GCoord> computeRotation(vector<GCoord>& v, XForm& M) {
  vector<GCoord> u;

  u.push_back(M.transform(v[0]));
  u.push_back(M.transform(v[1]));
  u.push_back(M.transform(v[2]));

  return(u);
}


double computeError(vector<GCoord> a, vector<GCoord> b) {
  double d=0;

  for (int i=0; i<3; i++) {
    GCoord c = matchSigns(a[i], b[i]);
    d += a[i].distance(c);
  }
  d /= 3.0;
  return(d);
}



int main() {
  int iter;
  int coarse_failures = 0;
  int fine_failures = 0;

  for (iter = 0; iter < maxcount; iter++) {
    AtomicGroup atoms = createGroup(nparticles, ellipsoid[0], ellipsoid[1], ellipsoid[2]);
    XForm M = randomXForm();

    vector<GCoord> directions = computeRotation(M);

    atoms.applyTransform(M);
    M.pop();
    atoms.applyTransform(M);
    vector<GCoord> axes = atoms.principalAxes();

    if (show_results) {
      cout << "Principal axes: [" <<axes[3] << "]\t" << axes[0] << "\t" << axes[1] << "\t" << axes[2] << endl;
      cout << " Computed axes: \t\t\t\t" << directions[0] << "\t" << directions[1] << "\t" << directions[2] << endl;
    }

    double d = computeError(axes, directions);
    if (show_results) {
      cout << "          ====> " << d << endl;
      cout << endl;
    }

    if (d >= threshold) {
      ++coarse_failures;
      cout << "***ERROR*** Failure (" << d << ") in self-check with threshold " << threshold << " at iteration " << iter << endl;
      if (exit_on_failure)
	exit(-1);
    }

    // ------------------------------------

    M = randomXForm();
    vector<GCoord> comp_dirs = computeRotation(axes, M);
    atoms.applyTransform(M);
    M.pop();
    atoms.applyTransform(M);
    vector<GCoord> new_axes = atoms.principalAxes();

    if (show_results) {
      cout << "Principal axes: [" <<new_axes[3] << "]\t" << new_axes[0] << "\t" << new_axes[1] << "\t" << new_axes[2] << endl;
      cout << " Original axes: \t\t\t\t" << comp_dirs[0] << "\t" << comp_dirs[1] << "\t" << comp_dirs[2] << endl;
    }

    d = computeError(new_axes, comp_dirs);
    if (show_results) {
      cout << "          ====> " << d << endl;
      cout << endl;
    }

    if (d >= fine_threshold) {
      ++fine_failures;
      cout << "***ERROR*** Failure (" << d << ") in self-check with fine-threshold " << fine_threshold << " at iteration " << iter << endl;
      if (exit_on_failure)
	exit(-1);
    }

  }

  cout << "There were " << setprecision(2) << coarse_failures*100.0/maxcount << "% failures in " << maxcount << " coarse tests.\n";
  cout << "There were " << setprecision(2) << fine_failures*100.0/maxcount << "% failures in " << maxcount << " fine tests.\n";

}
