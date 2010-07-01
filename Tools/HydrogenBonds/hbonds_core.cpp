/*
  Core code for hbonds tools

  (c) 2010 Tod D. Romo, Grossfield Lab, URMC
*/


#include <boost/format.hpp>

#include "hbonds_core.hpp"


double SimpleAtom::inner = 0.0;
double SimpleAtom::outer = 3.5;
double SimpleAtom::deviation = 20.0;



double SimpleAtom::distance2(const SimpleAtom& s) const
{
  double d;
  
  if (usePeriodicity)
    d = atom->coords().distance2(s.atom->coords(), sbox.box());
  else
    d = atom->coords().distance2(s.atom->coords());
  
  return(d);
}


double SimpleAtom::angle(const SimpleAtom& s) const {
  loos::GCoord left, middle, right;
  
  if (isHydrogen) {
    if (s.isHydrogen)
      throw(std::runtime_error("Cannot take the angle between two hydrogens"));
    
    left = attached_to->coords();
    middle = atom->coords();
    right = s.atom->coords();
    
  } else {
    if (!s.isHydrogen)
      throw(std::runtime_error("Cannot take the angle between two non-hydrogens"));
    
    left = atom->coords();
    middle = s.atom->coords();
    right = s.attached_to->coords();
  }
  
  if (usePeriodicity) {
    loos::GCoord box = sbox.box();
    
    left.reimage(box);
    middle.reimage(box);
    right.reimage(box);
  }
  
  return(loos::Math::angle(left, middle, right));
}



// Convert an AtomicGroup to a vector of SimpleAtoms.  We need the
// parent group (i.e. the system) in order to find atoms bound to
// atoms in the model group


std::vector<SimpleAtom> SimpleAtom::processSelection(const std::string& selection, const loos::AtomicGroup& system, const bool use_periodicity) {
  std::vector<SimpleAtom> processed;

  // We don't want to force model to be sorted (esp since it's const)
  // So, create a lightweight copy but use the same pAtoms...
  loos::AtomicGroup searchable(system);

  loos::AtomicGroup model = selectAtoms(system, selection);

  for (loos::AtomicGroup::const_iterator i= model.begin(); i != model.end(); ++i) {
    SimpleAtom new_atom(*i, system.sharedPeriodicBox(), use_periodicity);
    std::string n = (*i)->name();
    if (n[0] == 'H') {
      new_atom.isHydrogen = true;
      std::vector<int> bond_list = (*i)->getBonds();
      if (bond_list.size() == 0)
        throw(ErrorWithAtom(*i, "Detected a hydrogen bond that has no connectivity"));
      if (bond_list.size() != 1)
        throw(ErrorWithAtom(*i, "Detected a hydrogen bond that has more than one atom bound"));

      loos::pAtom pa = searchable.findById(bond_list[0]);
      if (pa == 0)
        throw(ErrorWithAtom(*i, "Cannot find atom hydrogen is bound too"));

      new_atom.attached_to = pa;
    }

    processed.push_back(new_atom);
  }

  return(processed);
}



bool SimpleAtom::hydrogenBond(const SAtom& o) const {
  double dist = distance2(o);
  double angl = angle(o);

  if (dist >= inner && dist <= outer && fmod(fabs(angl - 180.0), 360.0) <= deviation)
    return(true);

  return(false);
}



loos::AtomicGroup SimpleAtom::findHydrogenBonds(const std::vector<SimpleAtom>& group, const bool findFirstOnly) {
  loos::AtomicGroup results;

  for (SAGroup::const_iterator i = group.begin(); i != group.end(); ++i) {
    bool b = hydrogenBond(*i);
    if (b) {
      results.append(i->atom);
      if (findFirstOnly)
        break;
    }
  }

  return(results);
}



BondMatrix SimpleAtom::findHydrogenBondsMatrix(const std::vector<SimpleAtom>& group, loos::pTraj& traj, loos::AtomicGroup& model, const uint maxt) const {

  if (maxt > traj->nframes()) {
    std::cerr << boost::format("Error- row clip (%d) exceeds trajectory size (%d)\n") % maxt % traj->nframes();
    exit(-10);
  }

  BondMatrix bonds(maxt, group.size());

  for (uint t = 0; t < maxt; ++t) {
    traj->readFrame(t);
    traj->updateGroupCoords(model);

    for (uint i = 0; i<group.size(); ++i)
      bonds(t, i) = hydrogenBond(group[i]);

  }

  return(bonds);
}

