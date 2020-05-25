import numpy as np
import itertools

from andreff.source import solver as AndreffSolver
from solver.ext_solver import solver
from liegroups import SE3, SO3
import utility.utils as utils

def eval_with_noise(my_solver, T_v_rel, T_c_rel, scale, num_tests, cons='RCH'):
    #Initialize Storage variables
    #Extrinsic calibration Errors
    dual_trans_error= np.empty((num_tests))
    rel_trans_error = np.empty((num_tests))
    dual_rot_error= np.empty((num_tests))
    rel_rot_error = np.empty((num_tests))

    #Scale errors
    dual_alpha_error = np.empty((num_tests))
    rel_alpha_error = np.empty((num_tests))

    #Optimization Results
    dual_primal = np.empty((num_tests))
    rel_primal = np.empty((num_tests))
    dual_gap = np.empty((num_tests))
    rel_gap = np.empty((num_tests))

    #Algorithm Evaluation
    dual_time = np.empty((num_tests))
    rel_time = np.empty((num_tests))

    #Sanity Check
    dual_gt_value = np.empty((num_tests))
    rel_gt_value = np.empty((num_tests))

    #Choose multiple subsets of the data
    for i in range(num_tests):
        print("Solving set: {}".format(i))

        one_result_dict = evaluate_results(my_solver, T_v_rel[i], T_c_rel[i], scale, cons)

        #Store the estimation errors
        print("Storing results")
        dual_time[i] = one_result_dict['dual_time']
        dual_gt_value[i] = one_result_dict['dual_gt_value']
        dual_primal[i] = one_result_dict['dual_primal']
        dual_gap[i] = one_result_dict['dual_gap']
        dual_trans_error[i] = one_result_dict['dual_trans_error']
        dual_rot_error[i] = one_result_dict['dual_rot_error']
        dual_alpha_error[i] = one_result_dict['dual_alpha_error']
        rel_time[i] = one_result_dict['rel_time']
        rel_gt_value[i] = one_result_dict['rel_gt_value']
        rel_primal[i] = one_result_dict['rel_primal']
        rel_gap[i] = one_result_dict['rel_gap']
        rel_trans_error[i] = one_result_dict['rel_trans_error']
        rel_rot_error[i] = one_result_dict['rel_rot_error']
        rel_alpha_error[i] = one_result_dict['rel_alpha_error']


    #Store data
    results_dict = {'dual_time':dual_time, 'dual_gt_value':dual_gt_value, 'dual_primal': dual_primal,
                    'dual_gap':dual_gap, 'dual_trans_error': dual_trans_error, 'dual_rot_error':dual_rot_error, 'dual_alpha_error': dual_alpha_error,
                    'rel_time':rel_time, 'rel_gt_value':rel_gt_value, 'rel_primal': rel_primal,
                    'rel_gap':rel_gap, 'rel_trans_error': rel_trans_error, 'rel_rot_error':rel_rot_error, 'rel_alpha_error': rel_alpha_error}

    return results_dict

def evaluate_results(my_solver, T_v_rel, T_c_rel, scale, cons):
    #Number to know that the solution failed
    fail = 1000

    #Extrinsic calibration inverse for evaluating error
    ext_cal_trans = SE3.from_matrix(my_solver.extcal).trans
    ext_cal_rot_inv = SE3.from_matrix(my_solver.extcal).rot.inv()

    # Set the pose data
    my_solver.set_T_v_rel(T_v_rel)
    my_solver.set_T_c_rel(T_c_rel, scale)

    # Solve for the extrinsic calibration
    dual_time, dual_gt_value, dual_primal, dual_gap, dual_opt, dual_solution, dualflag = my_solver.dual_solve(cons)
    rel_time, rel_gt_value, rel_primal, rel_gap, relax_opt, relax_solution, relflag = my_solver.relax_solve(cons)
    #and_extcal, and_scale, and_runtime =  my_solver.andreff_solver.solve()

    if dualflag: #Solution was found
        # Store the estimation errors
        dual_trans_error = np.linalg.norm(ext_cal_trans - dual_solution[:3, 3])
        dual_rot_error = np.linalg.norm(
            (ext_cal_rot_inv.dot(SO3.from_matrix(dual_solution[:3, :3], normalize=True))).log())
        dual_alpha_error = np.abs(scale - 1 / dual_opt[3])
        dual_results_dict = {'dual_time':dual_time, 'dual_gt_value':np.asscalar(dual_gt_value), 'dual_primal': np.asscalar(dual_primal),
                    'dual_gap':np.asscalar(dual_gap), 'dual_trans_error': dual_trans_error, 'dual_rot_error':dual_rot_error, 'dual_alpha_error': np.asscalar(dual_alpha_error)}
    else: #Solution was not found
        dual_results_dict = {'dual_time': fail, 'dual_gt_value': fail,
                             'dual_primal': fail,
                             'dual_gap': fail, 'dual_trans_error': fail,
                             'dual_rot_error': fail, 'dual_alpha_error': fail}
    if relflag:
        rel_trans_error = np.linalg.norm(ext_cal_trans - relax_solution[:3, 3])
        rel_rot_error = np.linalg.norm(
            (ext_cal_rot_inv.dot(SO3.from_matrix(relax_solution[:3, :3], normalize=True))).log())
        rel_alpha_error = np.abs(scale - 1 / relax_opt[3])
        rel_results_dict = {'rel_time': rel_time, 'rel_gt_value': np.asscalar(rel_gt_value),
                             'rel_primal': np.asscalar(rel_primal),
                             'rel_gap': np.asscalar(rel_gap), 'rel_trans_error': rel_trans_error,
                             'rel_rot_error': rel_rot_error, 'rel_alpha_error': np.asscalar(rel_alpha_error)}
    else: #Solution was not found
        rel_results_dict = {'rel_time': fail, 'rel_gt_value': fail,
                             'rel_primal': fail,
                             'rel_gap': fail, 'rel_trans_error': fail,
                             'rel_rot_error': fail, 'rel_alpha_error': fail}
    rel_results_dict.update(dual_results_dict)

    return rel_results_dict


def limit_constraints(my_solver, T_v_rel, T_c_rel, scale, num_cases, trans_sigma, rot_sigma, set_of_cons):
    #Ones that are required for the solution to work atm
    #set_of_cons = ['RH', 'RCH']
    constraint_results={}
    if trans_sigma > 0 or rot_sigma > 0:
        # Add noise to the measurements
        print("Adding noise to data")
        T_v_rel_noisy = [utils.add_noise(T_v_rel, trans_sigma, rot_sigma) for i in range(num_cases)]
        T_c_rel_noisy = [utils.add_noise(T_c_rel, trans_sigma, rot_sigma) for i in range(num_cases)]
        for cons in set_of_cons:
            constraint_results[cons] = eval_with_noise(my_solver, T_v_rel_noisy, T_c_rel_noisy, scale, num_cases, cons)
    else:
        for cons in set_of_cons:
            constraint_results[cons] = evaluate_results(my_solver, T_v_rel, T_c_rel, scale, cons)
    return constraint_results

def main():
    #Make all numpy errors stop the program
    np.seterr(all='raise')

    #Set important values
    offset = 14
    num_cases = 100
    brookshire = False

    if not brookshire:
        # Name the relevant filenames
        data_filename = "unscaled_pose_data.mat"
        #Load the data
        T_vki, T_ci, extcal =utils.load_data_1(data_filename)
        # Convert inertial poses to relative
        T_v_rel = utils.inertial_to_relative(T_vki, offset)
        T_c_rel = utils.inertial_to_relative(T_ci, offset)
        trans_sigma = 0.01
        rot_sigma = 0
        scale = 2
        data_filename = "./results/per{}t{}r_".format(int(trans_sigma * 1000), int(rot_sigma * 1000)) + data_filename
        set_of_cons = ['R', 'RC', 'RH', 'RCH']
        #set_of_cons = ['RCH']
    else:
        # Name the relevant filenames
        data_filename = "brookshire_data.mat"
        T_v_rel, T_c_rel, extcal = utils.load_brookshire_data(data_filename)
        scale = 1
        trans_sigma = 0
        rot_sigma = 0
        data_filename = "./results/brookshire_data.mat"
        set_of_cons = ['RCH']

    # Initialize solver
    my_solver = solver()
    my_solver.set_extcal(extcal)
    my_solver.set_scale(scale)

    results = limit_constraints(my_solver, T_v_rel, T_c_rel, scale, num_cases, trans_sigma, rot_sigma, set_of_cons)

    #Save the results
    utils.results_saver(data_filename, results)


    # #Choose a subset of the data
    # num_poses_total = len(T_v_rel)
    # num_poses = 10
    # subset_indices = np.random.choice(num_poses_total, num_poses, replace=False).tolist()
    # T_v_rel_sub = [T_v_rel[index] for index in subset_indices]
    # T_c_rel_sub = [T_c_rel[index] for index in subset_indices]
    # my_solver.set_T_v_rel(T_v_rel_sub)
    # my_solver.set_T_c_rel(T_c_rel_sub, 2)

    return

if __name__ == "__main__":
    main()
