/*
  trajinfo

  trajinfo [options] model trajectory
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009 Tod D. Romo
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
#include <boost/format.hpp>
#include <boost/tuple/tuple.hpp>

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;

string model_name, traj_name;
bool brief = false;
bool box_info = false;

string centroid_selection;



// @cond TOOLS_INTERNAL

class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() : brief(false), box_info(false), centroid_selection("") { }
  
  void addGeneric(opts::po::options_description& o) {
    o.add_options()
      ("brief", opts::po::value<bool>(&brief)->default_value(brief), "Minimal output")
      ("centroid", opts::po::value<string>(&centroid_selection), "Report average centroid")
      ("box", opts::po::value<bool>(&box_info)->default_value(box_info), "Report periodic box info");
  }

  string print() const {
    ostringstream oss;

    oss << boost::format("brief=%d,centroid='%s',box=%d") % brief % centroid_selection % box_info;
    return(oss.str());
  }

  bool brief, box_info;
  string centroid_selection;
};

// @endcond

typedef boost::tuple<GCoord, GCoord, GCoord, GCoord, GCoord> BoxInfo;

BoxInfo scanBoxes(pTraj& traj) {
  GCoord min, max, avg;
  GCoord mine, maxe;

  traj->rewind();
  traj->readFrame();
  mine = maxe = min = max = avg = traj->periodicBox();
  double minsize = min[0] * min[1] * min[2];
  double maxsize = max[0] * max[1] * max[2];

  while (traj->readFrame()) {
    GCoord box = traj->periodicBox();
    avg += box;
    double size = box[0] * box[1] * box[2];
    if (size < minsize) {
      minsize = size;
      min = box;
    }
    if (size > maxsize) {
      maxsize = size;
      max = box;
    }
    for (int i=0; i<3; ++i) {
      if (mine[i] > box[i])
        mine[i] = box[i];
      if (maxe[i] < box[i])
        maxe[i] = box[i];
    }
  }

  avg /= traj->nframes();
  BoxInfo result(avg, min, max, mine, maxe);
  return(result);
}



boost::tuple<GCoord, GCoord> scanCentroid(AtomicGroup& model, pTraj& traj) {
  vector<GCoord> centers;
  GCoord avg;

  traj->rewind();
  while (traj->readFrame()) {
    traj->updateGroupCoords(model);
    GCoord c = model.centroid();
    centers.push_back(c);
    avg += c;
  }

  avg /= traj->nframes();

  GCoord std;
  for (uint i=0; i<traj->nframes(); ++i) {
    GCoord c = centers[i] - avg;
    std += c*c;
  }

  std /= (traj->nframes() - 1);
  for (uint i=0; i<3; i++)
    std[i] = sqrt(std[i]);

  boost::tuple<GCoord, GCoord> result(avg, std);
  return(result);
}


uint verifyFrames(pTraj& traj) {
  uint n = 0;
  
  traj->rewind();
  while (traj->readFrame())
    ++n;

  return(n);
}


// Boost doesn't currently allow you to specify the width in an arg, so we have
// to build the format string to fake it...
const string fldpre("%20s: ");

void verbInfo(AtomicGroup& model, pTraj& traj, AtomicGroup& center, const bool centroid = false) {
  cout << boost::format(fldpre + "%s\n") % "Model name" % model_name;
  cout << boost::format(fldpre + "%s\n") % "Trajectory name" % traj_name;
  cout << boost::format(fldpre + "%d\n") % "Number of atoms" % traj->natoms();
  cout << boost::format(fldpre + "%d\n") % "Number of frames" % traj->nframes();
  uint n = verifyFrames(traj);
  cout << boost::format(fldpre + "%d\n") % "Actual frames" % n;

  cout << boost::format(fldpre + "%f\n") % "Timestep" % traj->timestep();
  if (traj->hasPeriodicBox()) {
    cout << boost::format(fldpre + "%s\n") % "Periodic box" % "yes";
    if (box_info) {
      BoxInfo box = scanBoxes(traj);
      cout << boost::format(fldpre + "%s\n") % "Average box" % boost::get<0>(box);
      cout << boost::format(fldpre + "%s\n") % "Smallest box" % boost::get<1>(box);
      cout << boost::format(fldpre + "%s\n") % "Largest box" % boost::get<2>(box);
      cout << boost::format(fldpre + "%s x %s\n") % "Box extents" % boost::get<3>(box) % boost::get<4>(box);
    }
  } else
    cout << boost::format(fldpre + "%s\n") % "Periodic box" % "no";

  if (centroid) {
    boost::tuple<GCoord, GCoord> res = scanCentroid(center, traj);
    cout << boost::format(fldpre + "%s +- %s\n") % "Average Centroid" % boost::get<0>(res) % boost::get<1>(res);
  }
}


void briefInfo(pTraj& traj) {
  cout << traj->natoms() << " " << traj->nframes() << " " << traj->timestep() << " " << traj->hasPeriodicBox() << endl;
}





int main(int argc, char *argv[]) {

  opts::BasicOptions* bopts = new opts::BasicOptions;
  opts::BasicTrajectory* tropts = new opts::BasicTrajectory;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(tropts).add(topts);

  if (!options.parse(argc, argv))
    exit(-1);

  if (tropts->skip != 0)
    cerr << "Warning:  --skip is ignored by this tool\n";

  AtomicGroup model = tropts->model;
  pTraj traj = tropts->trajectory;

  AtomicGroup center;
  if (!(topts->centroid_selection.empty()))
    center = selectAtoms(model, topts->centroid_selection);

  if (!brief)
    verbInfo(model, traj, center, !(topts->centroid_selection.empty()));
  else
    briefInfo(traj);
}

