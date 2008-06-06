/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
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



#include <ios>
#include <sstream>
#include <iomanip>

#include <assert.h>
#include <algorithm>

#include <boost/random.hpp>

#include <AtomicGroup.hpp>



AtomicGroup* AtomicGroup::clone(void) const {
  return(new AtomicGroup(*this));
}


AtomicGroup AtomicGroup::copy(void) const {
  ConstAtomIterator i;
  AtomicGroup res;

  for (i = atoms.begin(); i != atoms.end(); i++) {
    pAtom pa(new Atom(**i));
    res.append(pa);
  }
  res._sorted = _sorted;
  res.box = box;

  return(res);
}


// Internal: verify the index into the atom array...
int AtomicGroup::rangeCheck(int i) const {
  if (i < 0)
    i = atoms.size() + i;
  if ((unsigned int)i >= atoms.size())
    throw(out_of_range("Bad index for an atom"));

  return(i);
}


// Returns the ith atom...  
pAtom AtomicGroup::getAtom(const int i) const {
  int j = rangeCheck(i);

  return(atoms[j]);
}

pAtom& AtomicGroup::operator[](const int i) {
  int j = rangeCheck(i);
  return(atoms[j]);
}
  
// For const objects...
const pAtom& AtomicGroup::operator[](const int i) const {
  int j = rangeCheck(i);
  return(atoms[j]);
}

// Internal: removes an atom from this group based on the address of the shared pointer...
void AtomicGroup::deleteAtom(pAtom pa) {
  vector<pAtom>::iterator iter;

  iter = find(atoms.begin(), atoms.end(), pa);
  if (iter == atoms.end())
    throw(runtime_error("Attempting to delete a non-existent atom"));

  atoms.erase(iter);
}


// Append each atom from the passed vector onto this group...
void AtomicGroup::append(vector<pAtom> pas) {
  vector<pAtom>::iterator i;

  for (i=pas.begin(); i != pas.end(); i++)
    atoms.push_back(*i);
}


// Append all atoms from the passed group onto this one
void AtomicGroup::append(const AtomicGroup& grp) {
  vector<pAtom>::const_iterator i;

  for (i=grp.atoms.begin(); i != grp.atoms.end(); i++)
    addAtom(*i);
}


// Remove all atoms in the passed vector
void AtomicGroup::remove(vector<pAtom> pas) {
  vector<pAtom>::iterator i;

  for (i=pas.begin(); i != pas.end(); i++)
    deleteAtom(*i);
}


// Removes all atoms contained in the passed group from this one...
void AtomicGroup::remove(const AtomicGroup& grp) {
  vector<pAtom>::const_iterator i;

  for (i=grp.atoms.begin(); i != grp.atoms.end(); i++)
    deleteAtom(*i);
}


// Concatenation operations...
AtomicGroup& AtomicGroup::operator+=(const AtomicGroup& rhs) {
  append(rhs);
  return(*this);
}

AtomicGroup& AtomicGroup::operator+=(const pAtom& rhs) {
  atoms.push_back(rhs);
  return(*this);
}

AtomicGroup AtomicGroup::operator+(const AtomicGroup& rhs) {
  AtomicGroup res(*this);
  res += rhs;
  return(res);
}

AtomicGroup AtomicGroup::operator+(const pAtom& rhs) {
  AtomicGroup res(*this);
  res += rhs;
  return(res);
}


// Internal: sort the atom array by atomid
void AtomicGroup::sort(void) {
  CmpById comp;

  if (!sorted())
    std::sort(atoms.begin(), atoms.end(), comp);

  sorted(true);
}



// Internal: calculates the start and stop iterators given offset and len args
// as in PERL's substr()...

boost::tuple<AtomicGroup::AtomIterator, AtomicGroup::AtomIterator> AtomicGroup::calcSubsetIterators(const int offset, const int len) {
  unsigned int a, b;

  if (offset < 0) {
    b = atoms.size() + offset + 1;
    a = (len == 0) ? 0 : b - len;
  } else {
    a = offset;
    b = (len == 0) ? atoms.size() : a + len;
  }

  if (b-a >= atoms.size())
    throw(range_error("Indices out of bounds for subsetting"));

  boost::tuple<AtomIterator, AtomIterator> res(atoms.begin() + a, atoms.begin() + b);

  return(res);
}



AtomicGroup AtomicGroup::subset(const int offset, const int len) {
  AtomicGroup res;

  boost::tuple<AtomIterator, AtomIterator> iters = calcSubsetIterators(offset, len);
  res.atoms.insert(res.atoms.begin(), boost::get<0>(iters), boost::get<1>(iters));
  
  return(res);
}


AtomicGroup AtomicGroup::excise(const int offset, const int len) {
  AtomicGroup res;

  boost::tuple<AtomIterator, AtomIterator> iters = calcSubsetIterators(offset, len);

  res.atoms.insert(res.atoms.begin(), boost::get<0>(iters), boost::get<1>(iters));
  atoms.erase(boost::get<0>(iters), boost::get<1>(iters));

  return(res);
}


// Not a bright algorithm...

AtomicGroup AtomicGroup::intersect(const AtomicGroup& grp) {
  AtomicGroup res;

  vector<pAtom>::iterator j;
  vector<pAtom>::const_iterator i;

  for (j=atoms.begin(); j != atoms.end(); j++)
    for (i=grp.atoms.begin(); i != grp.atoms.end(); i++)
      if (*i == *j) {
	res.addAtom(*j);
	break;
      }
    
  return(res);
}


// Select atoms from the current group, adding them to a new group
// based on the sel functor/predicate...

AtomicGroup AtomicGroup::select(const AtomSelector& sel) {
  AtomicGroup res;

  vector<pAtom>::const_iterator i;
  for (i=atoms.begin(); i != atoms.end(); i++)
    if (sel(*i))
      res.addAtom(*i);

  return(res);
}


// Split up a group into a vector of groups based on unique segids...
vector<AtomicGroup> AtomicGroup::splitByUniqueSegid(void) const {
  ConstAtomIterator i;
  UniqueStrings unique;

  for (i = atoms.begin(); i != atoms.end(); i++)
    unique.add((*i)->segid());

  int n = unique.size();
  int j;
  vector<AtomicGroup> results(n);
  for (i = atoms.begin(); i != atoms.end(); i++) {
    j = unique.find((*i)->segid());
    if (j < 0)
      throw(runtime_error("Could not find an atom we already found..."));
	
    results[j].append(*i);
  }

  return(results);
}



// Find an atom based on atomid
// Returns 0 (null shared_ptr) if not found...
pAtom AtomicGroup::findById(const int id) {
  sort();
  int bottom = 0, top = size()-1, middle;

  while (top > bottom) {
    middle = bottom + (top - bottom) / 2;
    if (atoms[middle]->id() < id)
      bottom = middle + 1;
    else
      top = middle;
  }

  if (atoms[bottom]->id() == id)
    return(atoms[bottom]);

  return(pAtom());
}

//! Note: when calling this, you'll want to make sure you use the 
//! outermost group (eg the psf or pdb you used to create things, rather than
//! using a subselection, unless you're sure the subsection contains these
//! atoms as well.  The main use of this routine is to create a group of atoms
//! bound to another atom.
AtomicGroup AtomicGroup::groupFromID(const vector<int> &id_list) {
    AtomicGroup result;
    for (unsigned int i=0; i<id_list.size(); i++) {
        pAtom pa = findById(id_list[i]);
        if (!pa) throw(out_of_range("Atom id doesn't exist"));
        result.addAtom(pa);
    }
    return(result);
}

// Get all atoms associated with the residue that contains the
// passed atom...  The returned atoms will not be in order.  If
// you want that, then explicitly sort the group.
AtomicGroup AtomicGroup::getResidue(pAtom res) {
  AtomIterator i;
  AtomicGroup result;

  i = find(atoms.begin(), atoms.end(), res);
  if (i == atoms.end())
    return(result);

  AtomIterator j = i;

  while (j >= atoms.begin()) {
    if ((*j)->resid() == res->resid() && (*j)->segid() == res->segid())
      result.addAtom(*j);
    else
      break;
    --j;
  }

  j = i+1;
  while (j < atoms.end()) {
    if ((*j)->resid() == res->resid() && (*j)->segid() == res->segid())
      result.addAtom(*j);
    else
      break;

    ++j;
  }
  
  return(result);
}


// renumber the atomids of the group...
void AtomicGroup::renumber(const int start, const int stride) {
  AtomIterator i;
  int id = start;

  for (i=atoms.begin(); i != atoms.end(); i++, id += stride)
    (*i)->id(id);
}


// Get the min and max atomid's...
int AtomicGroup::minId(void) const {
  ConstAtomIterator i;

  if (atoms.size() == 0)
    return(-1);

  int min = atoms[0]->id();
  for (i = atoms.begin()+1; i != atoms.end(); i++)
    if ((*i)->id() < min)
      min = (*i)->id();

  return(min);
}


int AtomicGroup::maxId(void) const {
  ConstAtomIterator i;
  
  if (atoms.size() == 0)
    return(-1);
  int max = atoms[0]->id();
  for (i=atoms.begin()+1; i != atoms.end(); i++)
    if ((*i)->id() > max)
      max = (*i)->id();

  return(max);
}


int AtomicGroup::minResid(void) const {
  ConstAtomIterator i;

  if (atoms.size() == 0)
    return(-1);

  int min = atoms[0]->resid();
  for (i = atoms.begin()+1; i != atoms.end(); i++)
    if ((*i)->resid() < min)
      min = (*i)->resid();

  return(min);
}

int AtomicGroup::maxResid(void) const {
  ConstAtomIterator i;

  if (atoms.size() == 0)
    return(-1);

  int max = atoms[0]->resid();
  for (i = atoms.begin()+1; i != atoms.end(); i++)
    if ((*i)->resid() < max)
      max = (*i)->resid();

  return(max);
}


// Count the number of higher structural elements...
int AtomicGroup::numberOfResidues(void) const {

  if (atoms.size() == 0)
    return(0);

  ConstAtomIterator i;
  int n = 1;
  int curr_resid = atoms[0]->resid();

  for (i=atoms.begin()+1; i !=atoms.end(); i++)
    if ((*i)->resid() != curr_resid) {
      ++n;
      curr_resid = (*i)->resid();
    }

  return(n);
}


int AtomicGroup::numberOfChains(void) const {

  if (atoms.size() == 0)
    return(0);

  ConstAtomIterator i;
  int n = 1;
  string curr_segid = atoms[0]->segid();

  for (i=atoms.begin()+1; i !=atoms.end(); i++)
    if ((*i)->segid() != curr_segid) {
      ++n;
      curr_segid = (*i)->segid();
    }

  return(n);
}


// Bounding box for all atoms in this group
vector<GCoord> AtomicGroup::boundingBox(void) const {
  greal min[3] = {0,0,0}, max[3] = {0,0,0};
  ConstAtomIterator i;
  int j;
  vector<GCoord> res(2);
  GCoord c;

  if (atoms.size() == 0) {
    res[0] = c;
    res[1] = c;
    return(res);
  }

  for (j=0; j<3; j++)
    min[j] = max[j] = (atoms[0]->coords())[j];

  for (i=atoms.begin()+1; i != atoms.end(); i++)
    for (j=0; j<3; j++) {
      if (max[j] < ((*i)->coords())[j])
	max[j] = ((*i)->coords())[j];
      if (min[j] > ((*i)->coords())[j])
	min[j] = ((*i)->coords())[j];
    }

  c.set(min[0], min[1], min[2]);
  res[0] = c;
  c.set(max[0], max[1], max[2]);
  res[1] = c;

  return(res);
}

// Geometric center of the group
GCoord AtomicGroup::centroid(void) const {
  GCoord c(0,0,0);
  ConstAtomIterator i;

  for (i = atoms.begin(); i != atoms.end(); i++)
    c += (*i)->coords();

  c /= atoms.size();
  return(c);
}


GCoord AtomicGroup::centerOfMass(void) const {
  GCoord c(0,0,0);
  ConstAtomIterator i;

  for (i=atoms.begin(); i != atoms.end(); i++) {
    c += (*i)->mass() * (*i)->coords();
  }
  c /= totalMass();
  return(c);
}


GCoord AtomicGroup::centerOfCharge(void) const {
  GCoord c(0,0,0);
  ConstAtomIterator i;

  for (i=atoms.begin(); i != atoms.end(); i++) {
    c += (*i)->charge() * (*i)->coords();
  }
  c /= totalCharge();
  return(c);
}

GCoord AtomicGroup::dipoleMoment(void) const {
    GCoord center = centerOfCharge();
    GCoord moment(0,0,0);
    ConstAtomIterator i;
    for (i=atoms.begin(); i != atoms.end(); i++) {
        moment += (*i)->charge() * ((*i)->coords() - center);
    }
    return(moment);
}

greal AtomicGroup::totalCharge(void) const {
  ConstAtomIterator i;
  greal charge = 0.0;

  for (i = atoms.begin(); i != atoms.end(); i++)
    charge += (*i)->charge();

  return(charge);
}

greal AtomicGroup::totalMass(void) const {
  ConstAtomIterator i;
  greal mass = 0.0;

  for (i = atoms.begin(); i != atoms.end(); i++)
    mass += (*i)->mass();

  return(mass);
}

// Geometric max radius of the group (relative to the centroid)
greal AtomicGroup::radius(void) const {
  GCoord c = centroid();
  greal radius = 0.0;
  ConstAtomIterator i;

  for (i=atoms.begin(); i != atoms.end(); i++) {
    greal d = c.distance2((*i)->coords());
    if (d > radius)
      radius = d;
  }

  radius = sqrt(radius);
  return(radius);
}



greal AtomicGroup::radiusOfGyration(void) const {
  GCoord c = centerOfMass();
  greal radius = 0;
  ConstAtomIterator i;

  for (i = atoms.begin(); i != atoms.end(); i++)
    radius += c.distance2((*i)->coords());

  radius = sqrt(radius / atoms.size());
  return(radius);
}


greal AtomicGroup::rmsd(AtomicGroup& v) {
  
  if (size() != v.size())
    throw(runtime_error("Cannot compute RMSD between groups with different sizes"));

  sort();
  v.sort();

  int n = size();
  double d = 0.0;
  for (int i = 0; i < n; i++) {
    GCoord x = atoms[i]->coords();
    GCoord y = v.atoms[i]->coords();
    d += x.distance2(y);
  }
  
  d = sqrt(d/n);

  return(d);
}



AtomicGroup operator+(const pAtom& lhs, const pAtom& rhs) {
  AtomicGroup res;
  res.append(lhs);
  res.append(rhs);
  return(res);
}


AtomicGroup operator+(const pAtom& lhs, const AtomicGroup& rhs) {
  AtomicGroup res(rhs);
  res += lhs;
  return(res);
}



bool AtomicGroup::operator==(AtomicGroup& rhs) {

  if (size() != rhs.size())
    return(false);

  if (this == &rhs)
    return(true);

  sort();
  rhs.sort();

  int n = size();
  for (int i = 0; i < n; i++)
    if (atoms[i] != rhs.atoms[i])
      return(false);

  return(true);
}


bool AtomicGroup::operator==(const AtomicGroup& rhs) const {
  if (size() != rhs.size())
    return(false);

  if (this == &rhs)
    return(true);

  const vector<pAtom> *lp;
  const vector<pAtom> *rp;
  vector<pAtom> lhs_atoms, rhs_atoms;
  CmpById comp;
  if (!sorted()) {
    lhs_atoms = atoms;
    std::sort(lhs_atoms.begin(), lhs_atoms.end(), comp);
  } else
    lp = &atoms;

  if (!rhs.sorted()) {
    rhs_atoms = rhs.atoms;
    std::sort(rhs_atoms.begin(), rhs_atoms.end(), comp);
    rp = &rhs_atoms;
  } else
    rp = &rhs.atoms;

  int n = size();
  for (int i = 0; i<n; i++)
    if ((*lp)[i] != (*rp)[i])
      return(false);

  return(true);
}


// XMLish output...
ostream& operator<<(ostream& os, const AtomicGroup& grp) {
  AtomicGroup::ConstAtomIterator i;
  if (grp._periodic)
    os << "<GROUP PERIODIC='" << grp.box << "'>\n";
  else
    os << "<GROUP>\n";
  for (i=grp.atoms.begin(); i != grp.atoms.end(); i++)
    os << "   " << **i << endl;
  os << "</GROUP>";

  return(os);
}



vector<GCoord> AtomicGroup::getTransformedCoords(const XForm& M) const {
  vector<GCoord> crds(atoms.size());
  ConstAtomIterator i;
  GMatrix W = M.current();
  int j = 0;

  for (i = atoms.begin(); i != atoms.end(); i++) {
    GCoord res = W * (*i)->coords();
    crds[j++] = res;
  }

  return(crds);
}



void AtomicGroup::applyTransform(const XForm& M) {
  AtomIterator i;
  GMatrix W = M.current();

  for (i = atoms.begin(); i != atoms.end(); i++)
    (*i)->coords() = W * (*i)->coords();

}



void AtomicGroup::dumpMatrix(const string s, double* A, int m, int n) const {
  cout << s << " = [\n" << endl;
  int i, j;

  for (j=0; j<m; j++) {
    for (i=0; i<n; i++)
      cout << setw(10) << A[i*m+j] << " ";
    cout << ";\n";
  }
  cout << "];\n";

}

// Returns a newly allocated array of double coords in row-major
// order...
double* AtomicGroup::coordsAsArray(void) const {
  double *A;
  int n = size();

  A = new double[n*3];
  int k = 0;
  int i;
  for (i=0; i<n; i++) {
    A[k++] = atoms[i]->coords().x();
    A[k++] = atoms[i]->coords().y();
    A[k++] = atoms[i]->coords().z();
  }

  return(A);
}

// Returns a newly allocated array of double coords in row-major order
// transformed by the current transformation.
double* AtomicGroup::transformedCoordsAsArray(const XForm& M) const {
  double *A;
  GCoord x;
  int n = size();
  GMatrix W = M.current();

  A = new double[n*3];
  int k = 0;
  int i;
  for (i=0; i<n; i++) {
    x = W * atoms[i]->coords();
    A[k++] = x.x();
    A[k++] = x.y();
    A[k++] = x.z();
  }

  return(A);
}



GCoord AtomicGroup::centerAtOrigin(void) {
  GCoord c = centroid();
  AtomIterator i;

  for (i = atoms.begin(); i != atoms.end(); i++)
    (*i)->coords() -= c;

  return(c);
}



void AtomicGroup::perturbCoords(const greal rms) {
  int i, n = size();
  GCoord r;

  loos::base_generator_type& rng = loos::rng_singleton();
  boost::uniform_real<> uni;
  boost::variate_generator<loos::base_generator_type&, boost::uniform_real<> > func(rng, uni);

  for (i=0; i<n; i++) {
    r.x(func());
    r.y(func());
    r.z(func());

    r /= r.length();
    r *= rms;

    atoms[i]->coords() += r;
  }
}


#if defined(__linux__) || defined(__APPLE__)


vector<GCoord> AtomicGroup::principalAxes(void) const {
  // Extract out the group's coordinates...
  int i;
  int n = size();
  double M[3] = {0.0, 0.0, 0.0};
  int k = 0;

  double *A = coordsAsArray();
  for (i=0; i<n; i++) {
    M[0] += A[k++];
    M[1] += A[k++];
    M[2] += A[k++];
  }

  M[0] /= n;
  M[1] /= n;
  M[2] /= n;

  // Subtract off the mean...
  for (i=k=0; i<n; i++) {
    A[k++] -= M[0];
    A[k++] -= M[1];
    A[k++] -= M[2];
  }

  // Multiply A*A'...
  double C[9];
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
	      3, 3, n, 1.0, A, 3, A, 3, 0.0, C, 3);

  delete A;

  // Now compute the eigen-decomp...
  char jobz = 'V', uplo = 'U';
  f77int nn;
  f77int lda = 3;
  double W[3], work[128];
  f77int lwork = 128;   // ???  Just a guess for sufficient storage to be
		     // efficient... 
  f77int info;
  nn = 3;

  dsyev_(&jobz, &uplo, &nn, C, &lda, W, work, &lwork, &info);
  if (info < 0)
    throw(runtime_error("dsyev_ reported an argument error..."));

  if (info > 0)
    throw(runtime_error("dsyev_ failed to converge..."));

  vector<GCoord> results(4);
  GCoord c;

  k = 0;
  for (i=0; i<3; i++) {
    c[0] = C[k++];
    c[1] = C[k++];
    c[2] = C[k++];
    results[2-i] = c;
  }

  // Now push the eigenvalues on as a GCoord...
  c[0] = W[2];
  c[1] = W[1];
  c[2] = W[0];

  results[3] = c;

  return(results);
}




GMatrix AtomicGroup::superposition(const AtomicGroup& grp) {
  int i, j;
  XForm W;

  int n = size();
  // Center both groups at the origin...

  GCoord xc = centroid();
  W.translate(-xc);
  double *X = transformedCoordsAsArray(W);
  

  GCoord yc = grp.centroid();
  W.identity();
  W.translate(-yc);
  double *Y = grp.transformedCoordsAsArray(W);

  // Compute correlation matrix...
  double R[9];
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 3, 3, n, 1.0, X, 3, Y, 3, 0.0, R, 3);

  double det = R[0]*R[4]*R[8] + R[3]*R[7]*R[2] + R[6]*R[1]*R[5] -
    R[0]*R[7]*R[5] - R[3]*R[1]*R[8] - R[6]*R[4]*R[2];

  // Now compute the SVD of R...
  char jobu = 'A', jobvt = 'A';
  f77int m = 3, lda = 3, ldu = 3, ldvt = 3, lwork=100, info;
  double work[lwork];
  f77int nn = 3;
  double S[3], U[9], Vt[9];
  
  dgesvd_(&jobu, &jobvt, &m, &nn, R, &lda, S, U, &ldu, Vt, &ldvt, work, &lwork, &info);

  // Adjust U (if necessary)
  if (det < 0.0) {
    U[6] = -U[6];
    U[7] = -U[7];
    U[8] = -U[8];
  }

  // Compute the rotation matrix...
  double M[9];
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, U, 3, Vt, 3, 0.0, M, 3);

  // Construct the new transformation matrix...  (W = M')
  GMatrix Z;
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      Z(i,j) = M[i*3+j];

  //W(0,3) = yc.x();
  //W(1,3) = yc.y();
  //W(2,3) = yc.z();


  W.identity();
  W.translate(yc);
  W.concat(Z);
  W.translate(-xc);

  delete[] X;
  delete[] Y;


  return(W.current());
}


GMatrix AtomicGroup::alignOnto(const AtomicGroup& grp) {
  XForm W;
  GMatrix M = superposition(grp);

  W.load(M);
  applyTransform(W);

  return(M);
}



#endif   /* defined(__linux__) || defined(__APPLE__) */
