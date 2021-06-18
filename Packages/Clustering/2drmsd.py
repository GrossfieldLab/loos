#!/usr/bin/env python3

"""
2drmsd.py: Calcuate 2D RMSD using different alignment and rmsd selections

rmsds aligns and calculates for same selection. The purpose of this tool is to calcuate 2D RMSD for a selection after aligning using a different selection. This script supports upto 2 trajectories.
"""



import loos
import argparse
import loos.pyloos
import math
import numpy

parser = argparse.ArgumentParser(description="Calcuate 2D RMSD using different alignment and rmsd selections")

parser.add_argument('model', help="Model")

parser.add_argument('traj', help="Trajectory")

parser.add_argument('--model2', help="Model for 2nd trajectory (default = model)" , default= None)

parser.add_argument('--traj2', help="2nd trajectory (default = traj)" , default= None)

parser.add_argument('--align', help="Align selection 1 for 1st trajectory (default = CA)", default= 'name == "CA"', type = str)

parser.add_argument('--align2', help="Align selection 2 for 2nd trajectory (default = align)", default=None, type = str)

parser.add_argument('--rmsd', help="RMSD selection for 1st trajectory (default = (C|O|N|CA))", default='name =~ "^(C|O|N|CA)$"', type = str)

parser.add_argument('--rmsd2', help="RMSD selection for 2nd trajectory (default = rmsd) ", default=None, type = str)

parser.add_argument('--precision', help="Write out matrix coefficients with this many digits (default = 2)", default=2, type = int)

parser.add_argument('--range1', help="Range of frames to use from the first trajectory [FORMAT - skip:stride:stop, stop is optional] (default = 0:1:last)", default='0:1', type = str)

parser.add_argument('--range2', help="Range of frames to use from the second trajectory [FORMAT - skip:stride:stop, stop is optional] (default = range1)", default=None,  type = str)

args=parser.parse_args()

#In case of only one trajectory is given as input

if args.rmsd2 is None:
	args.rmsd2 = args.rmsd

if args.align2 is None:
	args.align2 = args.align

if args.model2 is None:
	args.model2=args.model

if args.traj2 is None:
	args.traj2=args.traj

if args.range2 is None:
	args.range2=args.range1

#Model and Trajectory

model = loos.createSystem(args.model)
traj = loos.createTrajectory(args.traj, model)
model2 = loos.createSystem(args.model2)
traj2 = loos.createTrajectory(args.traj2, model2)


#Organizing range of frames to consider in trajectories

range_1=args.range1.split(":")
range_2=args.range2.split(":")

# In case only skip and stride is given as input
if len(range_1)==2:
	skip_1=int(range_1[0])
	stride_1=int(range_1[1])
	stop_1=len(traj)
# In case only skip, stride and stop are given
elif len(range_1)==3:
	skip_1=int(range_1[0])
	stride_1=int(range_1[1])
	stop_1=int(range_1[2])
else:
	print("Error in range")
	exit(0)


# In case only skip and stride is given as input
if len(range_2)==2:
	skip_2=int(range_2[0])
	stride_2=int(range_2[1])
	stop_2=len(traj2)
# In case only skip, stride and stop are given
elif len(range_2)==3:
	skip_2=int(range_2[0])
	stride_2=int(range_2[1])
	stop_2=int(range_2[2])
else:
	print("Error in range")
	exit(0)



#Align and RMSD selection definitions
align_selection_1 = loos.selectAtoms(model, args.align)
align_selection_2 = loos.selectAtoms(model2, args.align2)

rmsd_selection_1 = loos.selectAtoms(model, args.rmsd)
rmsd_selection_2 = loos.selectAtoms(model2, args.rmsd2)


print("# Align1 - ", args.align)
print("# Align2 - ", args.align2)
print("# RMSD1 - ", args.rmsd)
print("# RMSD2 - ", args.rmsd2)
print("# Traj-range1 - ", skip_1, ":" , stride_1, ":", stop_1)
print("# Traj-range2 - ", skip_2, ":" , stride_2, ":", stop_2)


for i in range(skip_1,stop_1,stride_1):
	traj.readFrame(i)
	traj.updateGroupCoords(model)
	ref_align=align_selection_1.copy()
	ref_target=rmsd_selection_1.copy()

	for j in range(skip_2,stop_2,stride_2):
		traj2.readFrame(j)
		traj2.updateGroupCoords(model2)
		# Find transformation matrix that aligns align_selection_2 onto ref_align(align_selection_1)
		trans_matrix = ref_align.alignOnto(align_selection_2)
		# Tranform matrix
		transform = loos.XForm(trans_matrix)
		# Apply the tranform to ref_target(rmsd_selection_1)
		ref_target.applyTransform(transform)
		# Calculate the RMSD between rmsd_selection_2 and ref_target(rmsd_selection_1)
		rmsd_value = rmsd_selection_2.rmsd(ref_target)
		# Print RMSD
		print(round(rmsd_value,args.precision),"",end='')

	print("")
