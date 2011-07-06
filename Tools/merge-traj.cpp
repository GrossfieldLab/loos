/*
  merge-dcd: combine multiple trajectories into a single long trajectory.  If the target
             trajectory exists, append to it.

  Alan Grossfield
  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School
 
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Alan Grossfield
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

using namespace std;

using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

// global for parsing program options
string model_name, output_traj, output_traj_downsample;
string center_selection;
vector<string> input_dcd_list;
int downsample_rate;
bool skip_first_frame=false;
bool reimage_by_molecule=false;


// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage
{
public:

  void addGeneric(po::options_description& o)
  {
    o.add_options()
      ("downsample-dcd", po::value<string>(&output_traj_downsample),
       "Downsampled DCD, must be synced with output_traj")
      ("downsample-rate", po::value<int>(&downsample_rate)->default_value(10),
       "Write every nth frame to downsampled DCD")
      ("centering-selection", po::value<string>(&center_selection)->default_value(""), "Selection for centering")
      ("skip-first-frame", po::value<bool>(&skip_first_frame)->default_value(false), "Skip first frame of each trajectory (for xtc files)")
      ("fix-imaging", po::value<bool>(&reimage_by_molecule)->default_value(false), "Reimage the system so molecules aren't broken across image boundaries")
      ;
  }

  string print() const {
    ostringstream oss;

    oss << boost::format("downsample-dcd='%s', downsample-rate=%d, centering-selection='%s', skip-first-frame=%d, fix-imaging=%d")
      % output_traj_downsample
      % downsample_rate
      % center_selection
      % skip_first_frame
      % reimage_by_molecule;

    return(oss.str());
  }
};

// @endcond




int main(int argc, char *argv[])
    {

    string hdr = invocationHeader(argc, argv);
    opts::BasicOptions* bopts = new opts::BasicOptions;
    ToolOptions* topts = new ToolOptions;
    opts::RequiredArguments* ropts = new opts::RequiredArguments;
    ropts->addArgument("model", "model-filename");
    ropts->addArgument("output_traj", "output-dcd");
    ropts->addVariableArguments("input_traj", "trajectory");

    opts::AggregateOptions options;
    options.add(bopts).add(topts).add(ropts);
    if (!options.parse(argc, argv))
      exit(-1);

    model_name = ropts->value("model");
    output_traj = ropts->value("output_traj");
    input_dcd_list = ropts->variableValues("input_traj");

    cout << hdr << endl;
    AtomicGroup system = createSystem(model_name);
    bool do_recenter = true;
    if ( center_selection.length() == 0 )
        {
        do_recenter = false;
        }

    DCDWriter output(output_traj, true);

    DCDWriter *output_downsample = 0;
    bool do_downsample = (output_traj_downsample.length() > 0);
    if (do_downsample)
        {
        output_downsample = new DCDWriter(output_traj_downsample, true);
        }

    // Set up to do the recentering
    vector<AtomicGroup> molecules;
    AtomicGroup center;
    vector<AtomicGroup>::iterator m;
    if ( do_recenter )
        {
        center = selectAtoms(system, center_selection);
        }

    if ( do_recenter || reimage_by_molecule)
        {
        if ( system.hasBonds() )
            {
            molecules = system.splitByMolecule();
            }
        else
            {
            molecules = system.splitByUniqueSegid();
            }
        }

    uint original_num_frames = output.framesWritten();
    cout << "Target trajectory " 
         << output_traj
         << " has " 
         << original_num_frames
         << " frames."
         << endl;

    uint previous_frames = 0;
    vector<string>::iterator f;
    for (f=input_dcd_list.begin(); f!=input_dcd_list.end(); ++f)
        {
        pTraj traj=createTrajectory(*f, system);
        int nframes = traj->nframes();
        if (skip_first_frame && nframes > 1)
            {
            nframes--;
            }
        cout << "File: " << *f << ": " << nframes;

        if ( previous_frames + nframes <= original_num_frames) 
            // all of this file is contained in the existing file, skip it
            {
            // increment the frame pointer
            previous_frames += nframes;
            cout << " ( " << previous_frames << " )"
                 << "\tSkipping trajectory " 
                 << endl;

            }
        else
            // we need at least some of the data from this file
            {
            int frames_to_skip = original_num_frames - previous_frames;
            if ( frames_to_skip > 0 )
                {
                traj->seekFrame(frames_to_skip-1);
                }
            else
                {
                frames_to_skip = 0;
                }

            // if this is an xtc file, we need to skip 1 more frame
            if (skip_first_frame)
                {
                traj->readFrame();
                }

            cout << " ( " << previous_frames + nframes - frames_to_skip
                 << " ) "
                 << "\t Writing " << nframes - frames_to_skip 
                 << " frames."
                 << endl;

            while ( traj->readFrame() )
                {
                traj->updateGroupCoords(system);

                // Find the smallest box dimension
                GCoord box = system.periodicBox();
                double smallest=1e20;
                for (int i=0; i<3; i++)
                    {
                    if (box[i] < smallest)
                        {
                        smallest = box[i];
                        }
                    }

                smallest /=2.0;


                // If molecules can be broken across image bondaries
                // (eg GROMACS), then we may need 2 translations to 
                // fix them -- first, translate the whole molecule such 
                // that a single atom is at the origin, reimage the
                // molecule, and put it back
                if (reimage_by_molecule)
                    {
                    for (m=molecules.begin(); m != molecules.end(); ++m )
                        {
                        // This is relatively slow, so we'll skip the 
                        // cases we know we won't need this -- 1 particle
                        // molecules and molecules with small radii 
                        if ( (m->size() > 1) && (m->radius() > smallest) )
                            {
                            m->mergeImage();
                            m->reimage();
                            }
                        }
                    }


                if ( do_recenter )
                    {
                    // Now, do the regular imaging.  Put the system centroid 
                    // at the origin, and reimage by molecule
                    GCoord centroid = center.centroid();
                    system.translate(-centroid);

                    for (m=molecules.begin(); m != molecules.end(); ++m )
                        {
                        m->reimage();
                        }

                    // Sometimes if the box has drifted enough, reimaging by molecule
                    // will significantly alter the centroid of the selected system, so
                    // we need to center a second time, which perversely means we'll need
                    // to reimage again. In my tests, this second go around is necessary and
                    // sufficient to fix everything, but I'm willing to be proved wrong.

                    centroid = center.centroid();
#if DEBUG
                    cerr << "centroid after reimaging: " << centroid << endl;
#endif
                    system.translate(-centroid);

                    for (m=molecules.begin(); m != molecules.end(); ++m )
                        {
                        m->reimage();
                        }
#if DEBUG
                    centroid = center.centroid();
                    cerr << "centroid after second reimaging: " << centroid << endl;
#endif 
                    }

                output.writeFrame(system);
                if ( do_downsample && (previous_frames % downsample_rate == 0) )
                    {
                    output_downsample->writeFrame(system);
                    }
                previous_frames++;
                }
            }

        }



    }
         
