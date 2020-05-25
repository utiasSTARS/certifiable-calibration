import numpy as np
import scipy.io as sio
from solver.ext_solver import solver
import utility.utils as utils

from liegroups import SE3

def gen_data():
    offset = 14
    num_poses = 50
    num_cases = 1

    # Name the relevant filenames
    data_filename = "unscaled_pose_data.mat"
    #Load the data
    T_vki, T_ci, extcal =utils.load_data_1(data_filename)
    # Convert inertial poses to relative
    T_v_rel = utils.inertial_to_relative(T_vki, offset)
    T_c_rel = utils.inertial_to_relative(T_ci, offset)
    trans_sigma = 0.1
    rot_sigma = 0.1
    scale = 2
    T_v_rel_noisy = utils.add_noise(T_v_rel, trans_sigma, rot_sigma)
    T_c_rel_noisy = utils.add_noise(T_c_rel, trans_sigma, rot_sigma)
    my_solver = solver()
    my_solver.set_T_c_rel(T_c_rel_noisy, scale)
    data_dict = {"T1": T_v_rel_noisy, "T2": my_solver.T_c_rel, "extcal": extcal, "scale":2}
    sio.savemat("failure_case.mat", data_dict)
    test = 1

def matlab_comparison():
    data = sio.loadmat('simple_case.mat')
    T1 = data['T1']
    T2 = data['T2']
    T_v_rel = [ SE3.from_matrix(T1[:, :, 0]), SE3.from_matrix(T1[:, :, 1])]
    T_c_rel = [ SE3.from_matrix(T2[:, :, 0]), SE3.from_matrix(T2[:, :, 1])]

    my_solver = solver()
    my_solver.set_T_v_rel(T_v_rel)
    my_solver.set_T_c_rel(T_c_rel)

    dual_time, dual_gt, dual_primal, dual_gap, dual_opt, dual_solution, dual_flag = my_solver.dual_solve('R', verbose=True)
    rel_time, rel_gt, rel_primal, rel_gap, relax_opt, relax_solution, rel_flag = my_solver.relax_solve('R')

    print(1/dual_opt[3])
    print(1/relax_opt[3])

    print(dual_solution)
    print(relax_solution)
    return

matlab_comparison()