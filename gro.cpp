#include <gro.hpp>



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
      float x = parseStringAs<float>(buf, 20, 8);
      float y = parseStringAs<float>(buf, 28, 8);
      float z = parseStringAs<float>(buf, 36, 8);
      
      // We ignore velocities...
      pAtom pa(new Atom);
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
      throw(std::runtime_error("Cannot parse box '" + buf + "'"));
    periodicBox(box);
  }
}
