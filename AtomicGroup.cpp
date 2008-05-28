/*
  AtomicGroup.cpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Basic class for groups of atoms...
*/


#include <ios>
#include <sstream>
#include <iomanip>

#include <assert.h>
#include <algorithm>

#include <AtomicGroup.hpp>


// Deep copy...  Creates new shared atoms and stuffs 'em
// into a new group...
AtomicGroup* AtomicGroup::clone(void) const {
  ConstAtomIterator i;
  AtomicGroup* res(new AtomicGroup);

  for (i = atoms.begin(); i != atoms.end(); i++) {
    pAtom pa(new Atom(**i));
    res->append(pa);
  }

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
      throw runtime_error("Could not find an atom we already found...");
	
    results[j].append(*i);
  }

  return(results);
}



// Need to check how find_if() actually searches since we have
// an ordered list...

AtomicGroup::AtomIterator AtomicGroup::findIteratorById(const int id) {
  AtomIterator i;

  if (!sorted())
    sort();

  BindId finder(id);
  i = find_if(atoms.begin(), atoms.end(), finder);
  return(i);
}


// Find an atom based on atomid
// Returns 0 (null shared_ptr) if not found...
pAtom AtomicGroup::findById(const int id) {
  AtomIterator i = findIteratorById(id);

  if (i == atoms.end())
    return(pAtom());

  return(*i);
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
AtomicGroup::BoundingBox AtomicGroup::boundingBox(void) const {
  BoundingBox bdd;
  ConstAtomIterator i;
  int j;

  if (atoms.size() == 0)
    return(bdd);

  for (j=0; j<3; j++)
    bdd.min[j] = bdd.max[j] = (atoms[0]->coords())[j];

  for (i=atoms.begin()+1; i != atoms.end(); i++)
    for (j=0; j<3; j++) {
      if (bdd.min[j] < ((*i)->coords())[j])
	bdd.min[j] = ((*i)->coords())[j];
      if (bdd.max[j] > ((*i)->coords())[j])
	bdd.max[j] = ((*i)->coords())[j];
    }

  return(bdd);
}

// Geometric center of the group
GCoord AtomicGroup::centroid(void) const {
  GCoord c(0,0,0);
  ConstAtomIterator i;

  for (i=atoms.begin(); i != atoms.end(); i++)
    c += (*i)->coords();

  c /= atoms.size();
  return(c);
}


GCoord AtomicGroup::centerOfMass(void) const {
  GCoord c(0,0,0);
  ConstAtomIterator i;

  for (i=atoms.begin(); i != atoms.end(); i++)
    c += (*i)->mass() * (*i)->coords();
  c /= atoms.size();
  return(c);
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


// XMLish output...
ostream& operator<<(ostream& os, const AtomicGroup& grp) {
  AtomicGroup::ConstAtomIterator i;
  os << "<GROUP>\n";
  for (i=grp.atoms.begin(); i != grp.atoms.end(); i++)
    os << "   " << **i << endl;
  os << "</GROUP>";

  return(os);
}



vector<GCoord> AtomicGroup::transformedCoords(void) const {
  vector<GCoord> crds(atoms.size());
  ConstAtomIterator i;
  int j = 0;

  for (i = atoms.begin(); i != atoms.end(); i++) {
    GCoord res = _xform.current() * (*i)->coords();
    crds[j++] = res;
  }

  return(crds);
}


void AtomicGroup::applyTransformation(void) {
  AtomIterator i;

  for (i = atoms.begin(); i != atoms.end(); i++) {
    GCoord res = _xform.current() * (*i)->coords();
    (*i)->coords(res);
  }
}



#if defined(__linux__) || defined(__APPLE__)

vector<GCoord> AtomicGroup::principalAxes(void) const {
  // Extract out the group's coordinates...
  int i;
  int n = size();
  double M[3] = {0.0, 0.0, 0.0};
  int k = 0;

  double *A = transformedCoordsAsArray();
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
  c[0] = W[0];
  c[1] = W[1];
  c[2] = W[2];

  results[3] = c;

  return(results);
}

#endif   /* defined(__linux__) || defined(__APPLE__) */


void AtomicGroup::dumpMatrix(const string s, double* A, int m, int n) const {
  cout << s << endl;
  int i, j;

  cout << "[" << endl;
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
double* AtomicGroup::transformedCoordsAsArray(void) const {
  double *A;
  GCoord x;
  int n = size();
  
  A = new double[n*3];
  int k = 0;
  int i;
  for (i=0; i<n; i++) {
    x = _xform.current() * atoms[i]->coords();
    A[k++] = x.x();
    A[k++] = x.y();
    A[k++] = x.z();
  }

  return(A);
}



void AtomicGroup::centerAtOrigin(void) {
  GCoord c = centroid();
  _xform.translate(-c);
}
