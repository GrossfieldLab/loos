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
#include <boost/program_options.hpp>

using namespace loos;
namespace po = boost::program_options;

// global for parsing program options
string model_name, output_traj, output_traj_downsample;
string center_selection;
vector<string> input_dcd_list;
int downsample_rate;
bool skip_first_frame;


void parseOptions(int argc, char *argv[])
    {
    try
        {
        po::options_description generic("Allowed options"); 
        generic.add_options()
            ("help", "Produce this message")
            ("downsample-dcd", po::value<string>(&output_traj_downsample),
                                    "Downsampled DCD, must be synced with output_traj")
            ("downsample-rate", po::value<int>(&downsample_rate)->default_value(10),
                                    "Write every nth frame to downsampled DCD")
            ("centering-selection", 
                    po::value<string>(&center_selection)->default_value(string("")),
                                    "Selection for centering")
            ("input_trajs", po::value< vector<string> >()->multitoken(), "Trajs to merge")
            ("skip-first-frame", "Skip first frame of each trajectory (for xtc files)")
            ;
            
        po::options_description hidden("Hidden options");
        hidden.add_options()
            ("model", po::value<string>(&model_name), "Model filename")
            ("output_traj", po::value<string>(&output_traj), "DCD filename")
            ;

        po::options_description command_line;
        command_line.add(hidden).add(generic);

        po::positional_options_description p;
        p.add("model", 1);
        p.add("output_traj", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc,argv).
                    options(command_line).positional(p).run(), vm);
        po::notify(vm);

        if ( vm.count("help") || !( vm.count("model") &&
                                    vm.count("output_traj") &&
                                    vm.count("input_trajs")
                                  )
           )                        
            {
            cerr << "Usage: " << argv[0]
                 << " model output_traj "
                 << "[options] "
                 << "--input_trajs file1 [file2 ...]"
                 << endl;
            cerr << generic;
            exit(-1);
            }

        input_dcd_list = vm["input_trajs"].as<vector<string> >();

        if (vm.count("skip-first-frame"))
            {
            skip_first_frame=true;
            }
        else
            {
            skip_first_frame=false;
            }


        }
    catch (exception &e)
        {
        cerr << "Error: " << e.what() << endl;
        exit(-1);
        }
    }


#if 0
void Usage()
    {
    cerr << "Usage: merge-dcd system recenter-selection output-dcdname "
            "input-dcd [input-dcd2 ...]"
         << endl;    
    cerr << "Giving a empty selection string turns off centering" << endl;
    cerr << "The input dcd files are concatenated in the command line order."
         << endl;
    }
#endif

int main(int argc, char *argv[])
    {
    parseOptions(argc, argv);
#if 0
    if ( (argc <=1) || 
         ( (argc >=2) && ((strncmp(argv[1], "-h", 2) == 0) ) ) ||
         (argc < 5) )
        {
        Usage();
        exit(-1);
        }
#endif

    cout << invocationHeader(argc, argv) << endl;
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

        if ( system.allHaveProperty(Atom::bondsbit) )
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
        if (skip_first_frame)
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

                if ( do_recenter )
                    {
                    GCoord centroid = center.centroid();
                    system.translate(-centroid);
                    for (m=molecules.begin(); m != molecules.end(); ++m )
                        {
                        m->reimage();
                        }
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
         
