// (c) 2012 Tod D. Romo, Grossfield Lab, URMC
// testing some ideas...

#include <iostream>
#include <string>
#include <cstdlib>

#include <netcdf.h>

using namespace std;

#if defined(__APPLE__)
typedef unsigned int     uint;
#endif


struct SetEdgesForCoords {
  static const uint ndims = 3;

  SetEdgesForCoords(const uint natoms) : _natoms(natoms) { }

  void setStart(size_t s[], const uint frame) const {
    s[0] = frame;
    s[1] = 0;
    s[2] = 0;
  }

  void setCount(size_t s[]) const {
    s[0] = 1;
    s[1] = _natoms;
    s[2] = 3;
  }

  size_t size() const { return(_natoms * 3); }


  uint _natoms;
};


struct SetEdgesForBoxes {
  static const uint ndims = 2;

  void setStart(size_t s[], const uint frame) const {
    s[0] = frame;
    s[1] = 0;
  }

  void setCount(size_t s[]) const {
    s[0] = 1;
    s[1] = 3;
  }

  size_t size() const { return(3); }
  
};



void dump(const uint n, const size_t s[], const string& msg) {
  cerr << msg << " = (";
  for (uint i=0; i<n; ++i)
    cerr << s[i] << (i == (n-1) ? ')' : ',');
  cerr << endl;
}

template<class SETTER>
struct Wrapper {
  Wrapper(const int ncid, const int varid, const SETTER& setter) :
    _ncid(ncid),
    _varid(varid),
    _setter(setter),
    _data(0)

  {
    // Get type...
    int retval = nc_inq_vartype(_ncid, _varid, &_type);

    // Allocate space
    switch(_type) {
    case NC_BYTE: _data = new unsigned char[_setter.size()]; break;
    case NC_CHAR: _data = new char[_setter.size()]; break;
    case NC_SHORT: _data = new short[_setter.size()]; break;
    case NC_INT: _data=new int[_setter.size()]; break;
    case NC_FLOAT: _data=new float[_setter.size()]; break;
    case NC_DOUBLE: _data=new double[_setter.size()]; break;
    }
  }


  template<typename T>
  T get(const uint i) {
    T x;

    switch(_type) {
    case NC_BYTE: x = (static_cast<unsigned char*>(_data))[i]; break;
    case NC_CHAR: x = (static_cast<char*>(_data))[i]; break;
    case NC_SHORT: x = (static_cast<short*>(_data))[i]; break;
    case NC_INT: x = (static_cast<int*>(_data))[i]; break;
    case NC_FLOAT: x = (static_cast<float*>(_data))[i]; break;
    case NC_DOUBLE: x = (static_cast<double*>(_data))[i]; break;
    }
    return(x);
  }
                        
  template<typename T>
  void set(const uint i, const T x) {
    switch(_type) {
    case NC_BYTE: x = (static_cast<unsigned char*>(_data))[i] = x; break;
    case NC_CHAR: x = (static_cast<char*>(_data))[i] = x; break;
    case NC_SHORT: x = (static_cast<short*>(_data))[i] = x; break;
    case NC_INT: (static_cast<int*>(_data))[i] = x; break;
    case NC_FLOAT: (static_cast<float*>(_data))[i] = x; break;
    case NC_DOUBLE: (static_cast<double*>(_data))[i] = x; break;
    }
  }

  nc_type type() const { return(_type); }


  bool readFrame(const uint frame) {
    _setter.setStart(_start, frame);
    _setter.setCount(_count);

    int retval = nc_get_vara(_ncid, _varid, _start, _count, _data);
    if (retval)
      cerr << "Internal error - " << retval << endl;
    return(retval == 0);
  }


  int _ncid, _varid;
  const SETTER _setter;
  void* _data;
  nc_type _type;
  size_t _start[SETTER::ndims];
  size_t _count[SETTER::ndims];
};








int main(int argc, char *argv[]) {
  int ncid;

  cout << "* Opening...\n";
  int retval;
  retval = nc_open(argv[1], NC_NOWRITE, &ncid);
  cout << "retval = " << retval << endl;
  cout << "ncid = " << ncid << endl;
  
  int ngatts, nvars;
  retval = nc_inq_natts(ncid, &ngatts);
  retval = nc_inq_nvars(ncid, &nvars);
  cout << "# global attribues: " << ngatts << endl;
  cout << "# variables: " << nvars << endl;

  cout << "* Global attributes:\n";
  for (uint i=0; i<ngatts; ++i) {
    char buf[512];
    retval = nc_inq_attname(ncid, NC_GLOBAL, i, buf);
    size_t l;
    retval = nc_inq_attlen(ncid, NC_GLOBAL, buf, &l);
    nc_type t;
    retval = nc_inq_atttype(ncid, NC_GLOBAL, buf, &t);
    cout << i << '\t' << buf << " :\t";
    if (t == NC_CHAR && l > 0) {
      char text[1024];
      retval = nc_get_att_text(ncid, NC_GLOBAL, buf, text);
      text[l] = '\0';
      cout << text;
    }
    cout << endl;
    
  }

  size_t badl;
  cout << "> Testing for missing attribute by length: ";
  retval = nc_inq_attlen(ncid, NC_GLOBAL, "snufkin", &badl);
  cout << retval << endl;
  


  int frame_id, spatial_id, atom_id, label_id, cell_spatial_id, cell_angular_id, cell_lengths_id;
  size_t len;

  retval = nc_inq_dimid(ncid, "frame", &frame_id);
  retval = nc_inq_dimlen(ncid, frame_id, &len);
  cout << "Frame len = " << len << endl;
  size_t nframes = len;

  retval = nc_inq_dimid(ncid, "atom", &atom_id);
  retval = nc_inq_dimlen(ncid, atom_id, &len);
  cout << "Atom len = " << len << endl;
  
  int coord_id;
  retval = nc_inq_varid(ncid, "coordinates", &coord_id);
  cout << "coord_id = " << coord_id << endl;
  int coord_ndims;
  retval = nc_inq_varndims(ncid, coord_id, &coord_ndims);
  cout << "coord_ndims = " << coord_ndims << endl;
  nc_type coord_type;
  retval = nc_inq_vartype(ncid, coord_id, &coord_type);
  cout << "Coord type = ";
  switch(coord_type) {
  case NC_FLOAT: cout << "float"; break;
  case NC_DOUBLE: cout << "double"; break;
  default: cout << "other";
  }
  cout << endl;

  int coord_natts;
  retval = nc_inq_varnatts(ncid, coord_id, &coord_natts);
  cout << "Coord Natts = " << coord_natts << endl;

  for (int i=0; i<coord_natts; ++i) {
    char buf[512];
    retval = nc_inq_attname(ncid, coord_id, i, buf);
    if (retval != 0)
      break;
    nc_type type;
    retval = nc_inq_atttype(ncid, coord_id, buf, &type);
    cout << '\t' << i << ": (";
    switch(type) {
    case NC_FLOAT: cout << "float"; break;
    case NC_DOUBLE: cout << "double"; break;
    case NC_CHAR: cout << "char"; break;
    case NC_INT: cout << "int"; break;
    default: cout << "other";
    }

    cout << ")\t" << buf << '\t';

    if (type == NC_CHAR) {
      char val[512];
      size_t l;
      retval = nc_inq_attlen(ncid, coord_id, buf, &l);
      retval = nc_get_att_text(ncid, coord_id, buf, val);
      val[l] = '\0';
      cout << "'" << val << "'";
    }
    cout << endl;
  }

  int coord_dimids[10];
  retval = nc_inq_vardimid(ncid, coord_id, coord_dimids);
  cout << "Dimids: ";
  for (int i=0; i<coord_ndims; ++i)
    cout <<  coord_dimids[i] << '\t';
  cout << endl;
  
  int coord_storage;
  size_t coord_chunksize;
  retval = nc_inq_var_chunking(ncid, coord_id, &coord_storage, &coord_chunksize);
  if (retval == 0) {
    cout << "storage = " << (coord_storage == NC_CONTIGUOUS ? "contig" : "chunked") << endl;
    if (coord_storage == NC_CHUNKED)
      cout << "chunk size = " << coord_chunksize << endl;
  } else
    cout << "No chunk info\n";


  cout << "*Reading first atom...\n";
  float first_atom[4] = {0,0,0,0};
  size_t first_start[3] = {0, 0, 0};
  size_t first_count[3] = {1, 1, 3};
  retval = nc_get_vara_float(ncid, coord_id, first_start, first_count, first_atom);
  cout << "retval = " << retval << endl;
  cout << "First coord = (" << first_atom[0] << ',' << first_atom[1] << ',' << first_atom[2] << ")\n";
  cout << "Probe = " << first_atom[3] << endl;

  cout << "*Reading a frame...\n";
  float *coords = new float[(len+1) * 3];
  for (size_t i =0; i<(len+1)*3; ++i)
    coords[i] = 0.0;

  size_t start[3] = {0,0,0};
  size_t count[3] = {1, 1, 3};
  count[1] = len;
  retval = nc_get_vara_float(ncid, coord_id, start, count, coords);
  cout << "retval = " << retval << endl;
  cout << "First atom = (" << coords[0] << ',' << coords[1] << ',' << coords[2] << ")\n";
  size_t idx = (len-1) * 3;
  cout << "Last atom = (" << coords[idx] << ',' << coords[idx+1] << ',' << coords[idx+2] << ")\n";
  cout << "Probe = " << coords[idx+3] << endl;

  cout << "\n*Reading 2nd frame...\n";
  for (size_t i =0; i<(len+1)*3; ++i)
    coords[i] = 0.0;

  start[0] = 1; start[1] = 0; start[2] = 0;
  count[0] = 1; count[1] = len; count[2] = 3;
  dump(3, start, "start");
  dump(3, count, "count");
  retval = nc_get_vara_float(ncid, coord_id, start, count, coords);
  cout << "retval = " << retval << endl;
  cout << "First atom = (" << coords[0] << ',' << coords[1] << ',' << coords[2] << ")\n";
  cout << "Last atom = (" << coords[idx] << ',' << coords[idx+1] << ',' << coords[idx+2] << ")\n";
  cout << "Probe = " << coords[idx+3] << endl;

  cout << "\n*Dump of periodic boxes...\n";
  retval = nc_inq_varid(ncid, "cell_lengths", &cell_lengths_id);
  if (retval) {
    cerr << "Failed.\n";
    exit(-1);
  }

  for (uint i=0; i<nframes; ++i) {
    size_t start[2]; start[0] = i; start[1] = 0;
    size_t count[2]; count[0] = 1; count[1] = 3;
    double box[3];
    retval = nc_get_vara_double(ncid, cell_lengths_id, start, count, box);
    cout << i << "\t(" << box[0] << ',' << box[1] << ',' << box[2] << ")\n";
  }

  
  cout << "\n* Wrapper interface test...\n";
  Wrapper<SetEdgesForCoords> wrc(ncid, coord_id, SetEdgesForCoords(len));
  wrc.readFrame(0);
  cout << "First atom = (" << wrc.get<float>(0) << ',' << wrc.get<float>(1) << ',' << wrc.get<float>(2) << ")\n";

}

