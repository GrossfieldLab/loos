/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo, Alan Grossfield
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



#include <gro.hpp>
#include <utils.hpp>
#include <Fmt.hpp>



namespace loos {
  
  void Gromacs::read(std::istream& ifs) {
    std::string buf;

    getline(ifs, title_);

    // Get the # of atoms;;;
    getline(ifs, buf);
    int natoms = parseStringAs<int>(buf);
    while (natoms-- > 0) {
      getline(ifs, buf);
      
      int resid = parseStringAs<int>(buf, 0, 5);
      std::string resname = parseStringAs<std::string>(buf, 5, 5);
      std::string name = parseStringAs<std::string>(buf, 10, 5);
      int atomid = parseStringAs<int>(buf, 15, 5);
      float x = parseStringAs<float>(buf, 20, 8) * 10.0;
      float y = parseStringAs<float>(buf, 28, 8) * 10.0;
      float z = parseStringAs<float>(buf, 36, 8) * 10.0;
      
      // We ignore velocities...
      pAtom pa(new Atom);
      pa->index(_max_index++);
      pa->resid(resid);
      pa->id(atomid);
      pa->resname(resname);
      pa->name(name);
      pa->coords(GCoord(x,y,z));

      append(pa);
    }

    // Now process box...
    getline(ifs, buf);
    std::istringstream iss(buf);
    GCoord box;
    if (!(iss >> box[0] >> box[1] >> box[2]))
      throw(FileReadError(_filename, "Cannot parse box '" + buf + "'"));
    periodicBox(box * 10.0);

    // Since the atomic field in .gro files is only 5-chars wide, it can
    // overflow.  if there are enough atoms to cause an overflow, manually
    // renumber everything...
    if (atoms.size() >= 100000)
      renumber();
  }

  std::string Gromacs::atomAsString(const pAtom p) const {
    std::ostringstream s;

    // Float formatter for coords
    Fmt crdfmt(3);
    crdfmt.width(8);
    crdfmt.right();
    crdfmt.trailingZeros(true);
    crdfmt.fixed();

    // Float formatter for velocities
    Fmt velfmt(4);
    velfmt.width(8);
    velfmt.right();
    velfmt.trailingZeros(true);
    velfmt.fixed();

    /*
    // Float formatter for B's and Q's...
    Fmt bqfmt(2);
    bqfmt.width(6);
    bqfmt.right();
    bqfmt.trailingZeros(true);
    bqfmt.fixed();

    // We don't worry about strings exceeding field-widths (yet),
    // but do check for numeric overflows...
    s << std::setw(6) << std::left << p->recordName();
    s << hybrid36AsString(p->id(), 5) << " ";
    s << std::setw(4) << std::left << p->name();

    s << std::setw(1) << p->altLoc();
    s << std::setw(4) << std::left << p->resname();
    
    s << std::setw(1) << std::right << p->chainId();
    s << hybrid36AsString(p->resid(), 4);
    s << std::setw(2) << p->iCode();
    s << "  ";
        
    s << crdfmt(p->coords().x());
    s << crdfmt(p->coords().y());
    s << crdfmt(p->coords().z());
    s << bqfmt(p->occupancy());
    s << bqfmt(p->bfactor());
    s << "      ";
    s << std::setw(4) << std::left << p->segid();
    s << std::setw(2) << std::right << p->PDBelement();
    if (_show_charge)
      s << std::setw(2) << p->charge();
    else
      s << "  ";
    */

    s << std::setw(5) << p->resid();
    s << std::left << std::setw(5) << p->resname();
    s << std::right << std::setw(5) << p->name();
    s << std::setw(5) << p->id();
    s << crdfmt(p->coords().x()/10.);
    s << crdfmt(p->coords().y()/10.);
    s << crdfmt(p->coords().z()/10.);
    s << velfmt(0.0);
    s << velfmt(0.0);
    s << velfmt(0.0);

    return(s.str());
    
  }

  //! Build a GRO from an AtomicGroup
  //*
  //  As it stands, this is a horrible hack and must be fixed.
  //  I did it this way because I couldn't figure out how to 
  //  make the following work: 
  //  Gromacs p(g);
  //
  Gromacs Gromacs::fromAtomicGroup(const AtomicGroup& g) {
    Gromacs p;
    for (uint i = 0; i < g.size(); ++i) {
        p.append(g[i]);
    }
    
    if (g.isPeriodic()) {
        p.periodicBox(g.periodicBox());
    }

    return(p);
  }

  

  //! Output the group as a GRO...

  std::ostream& operator<<(std::ostream& os, const Gromacs& g) {
      AtomicGroup::const_iterator i;
      
      os << g.title() << std::endl;
      os << g.size() << std::endl;
      for (i=g.atoms.begin(); i != g.atoms.end(); ++i) {
          os << g.atomAsString(*i) << std::endl;
      }

      GCoord box = g.periodicBox();
      box /= 10.0;
      os << box.x() << "  " << box.y() << "  " << box.z() << std::endl;
      
          
      return(os);
  }
}
