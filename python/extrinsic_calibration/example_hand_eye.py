import numpy as np
from solver.ext_solver import solver
from utility.utils import load_data_1, inertial_to_relative, add_noise

# Name the relevant filenames
data_filename = "unscaled_pose_data.mat"

# Load the data
T_vki, T_cki, extcal = load_data_1(data_filename)
scale = 2

#Set noise levels
trans_noise = 0.01 #Percentage of translation
rot_noise = 0.01 #Percentage of rotation

# Initialize solver
my_solver = solver()

# Convert inertial poses to relative
T_v_rel = inertial_to_relative(T_vki)
T_c_rel = inertial_to_relative(T_cki)

# Choose a subset of the data
num_poses_total = len(T_v_rel)
num_poses = 100
subset_indices = np.random.choice(num_poses_total, num_poses, replace=False).tolist()
T_v_rel_sub = add_noise([T_v_rel[index] for index in subset_indices], trans_noise, rot_noise)
T_c_rel_sub = add_noise([T_c_rel[index] for index in subset_indices], trans_noise, rot_noise)

#Load the data into the solver
my_solver.set_T_v_rel(T_v_rel_sub) #Load relative egomotion sensor poses
my_solver.set_T_c_rel(T_c_rel_sub, scale=scale) # Load and scale relative camera sensor poses

# Run Solver
dual_time, dual_gt, dual_primal, dual_gap, dual_opt, dual_solution, dual_flag = my_solver.dual_solve(cons="RCH")
rel_time, rel_gt, rel_primal, rel_gap, relax_opt, relax_solution, rel_flag = my_solver.relax_solve(cons="RCH")

print("Scale: {}".format(scale))
print("Estimated Scale: {}".format(1/dual_opt[3]))
print("Extrinsic Calibration:")
print(extcal)
print("Estimated Extrinsic Calibration:")
print(dual_solution)

test=1