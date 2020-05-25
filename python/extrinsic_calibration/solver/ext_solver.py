import time

import numpy as np
import solver.constraint_gen as cg
import cvxpy as cp
from utility.utils import load_data_1, inertial_to_relative

from andreff.source import solver as AndreffSolver

class solver():
    def __init__(self):
        self.T_v_rel=[]
        self.T_c_rel = []
        self.extcal = []
        self.gt = np.empty((14, 1))
        self.gt[-1, 0] = 1

        # Setup Andreff solver
        self.andreff_solver = AndreffSolver()

        return

    def set_T_v_rel(self, pose_list):
        self.T_v_rel = pose_list
        return

    def set_T_c_rel(self, pose_list, scale=1):
        self.T_c_rel = pose_list
        if scale is not 1:
            for i in range(len(self.T_c_rel)):
                self.T_c_rel[i].trans = scale * self.T_c_rel[i].trans
        return self.T_c_rel

    def set_AndreffSolver(self):
        self.andreff_solver.set_vals(self.T_v_rel, self.T_c_rel, self.extcal)
        return

    def set_scale(self, scale):
        '''Set this for error evaluation'''
        self.gt[12, 0] = scale
        return

    def set_extcal(self, extcal):
        ''' Set this for error evaluation'''
        self.extcal = extcal
        self.gt[:3, 0] = extcal[:3 ,3]
        self.gt[3:12, 0] = extcal[:3 ,:3].flatten('F')
        return

    def get_rot_matrix(self, pose_list):
        rot_matrices = [pose.rot.as_matrix() for pose in pose_list]
        return rot_matrices

    def get_trans_vec(self, pose_list):
        trans_list = [np.transpose(np.atleast_2d(pose.trans)) for pose in pose_list]
        return trans_list

    def construct_M_R(self):
        rot_a = self.get_rot_matrix(self.T_v_rel)
        rot_b = self.get_rot_matrix(self.T_c_rel)
        M_R_list = []
        for i in range(len(rot_a)):
            kron_prod = np.kron(np.transpose(rot_a[i]), np.eye(3)) - np.kron(np.transpose(np.eye(3)), rot_b[i])
            M_R_list.append(np.hstack((np.zeros((9, 4)), kron_prod, np.zeros((9, 1)))))
        return M_R_list

    def construct_M_t(self):
        M_t_list = []
        rot_a = self.get_rot_matrix(self.T_v_rel)
        rot_b = self.get_rot_matrix(self.T_c_rel)
        trans_a = self.get_trans_vec(self.T_v_rel)
        trans_b = self.get_trans_vec(self.T_c_rel)
        for i in range(len(rot_a)):
            term_1 = np.eye(3) - rot_b[i]
            term_2 = np.kron(np.transpose(trans_a[i]), np.eye(3))
            M_t_list.append(np.hstack((term_1, -trans_b[i], term_2, np.zeros((3, 1)))))
        return M_t_list

    def construct_Q_i(self, M_i_list):
        Q_i_i = np.empty((len(M_i_list), M_i_list[0].shape[1], M_i_list[0].shape[1]))
        for i, entry in enumerate(M_i_list):
            Q_i_i[i, :, :] = np.matmul(np.transpose(entry),entry)
        Q_i = np.sum(Q_i_i, axis=0)
        return Q_i

    def construct_Q(self, verbose=False):
        start = time.time()
        # Construct the rotation error cost
        M_R_list = self.construct_M_R()
        time_M_R = time.time() - start

        Q_R = self.construct_Q_i(M_R_list)
        time_Q_R = time.time() - start - time_M_R

        #Construct the translation error cost
        M_t_list = self.construct_M_t()
        time_M_t = time.time() - start - time_Q_R

        Q_t = self.construct_Q_i(M_t_list)
        time_Q_t = time.time() - start - time_M_t

        Q = Q_t + Q_R
        Q = (np.transpose(Q) + Q)/2 #Make sure that Q is symmetric

        if verbose:
            print('[time profile] Q generation: M_r {} secs, Q_r {} secs, M_t {} secs, Q_t {} secs'.format(time_M_R, time_Q_R, time_M_t, time_Q_t))
        return Q

    def dual_solve(self, cons='RCH', verbose=True):
        if verbose:
            print("Dual Relaxation")
        #Construct the weight matrix
        Q = self.construct_Q(verbose)
        if verbose:
            print(Q)

        #Check Rank of Q_tt
        if np.linalg.matrix_rank(Q[:3, :3]) < 3:
            print('Q_{t,t} rank < 3 for data given! Exiting Solver.')
            return None

        Q_tilde = Q[4:, 4:] - np.matmul(np.matmul(Q[4:, :4], np.linalg.inv(Q[:4, :4])), Q[:4, 4:])
        P, X_4 = cg.construct_P(constraints=cons)
        Z = Q_tilde + P
        constraint = [cp.PSD(Z)]
        obj = cp.Maximize(X_4)
        problem = cp.Problem(obj, constraint)

        #Solve the dual problem
        #problem.solve(verbose=True, solver=cp.SCS, eps=10**(-10))
        #Initialize timer
        start_time = time.time()
        problem.solve(verbose=verbose, solver=cp.CVXOPT, abstol=10**(-7), refinement=1, kktsolver='robust')
        runtime = time.time() - start_time

        #Exit if solution was not found
        if problem.status in ["infeasible", "unbounded"]:
            return runtime, 0, 0, 0, 0, 0, False

        #Calculate the dual value and the matrx Z
        dual_value = problem.value
        Z_opt = Z.value

        #Check Rank of Z_opt
        #Z_opt_rank = np.linalg.matrix_rank(Z_opt)
        U, S, _ = np.linalg.svd(Z_opt)
        nullspace = np.logical_or(np.isclose(np.log(S), np.log(S[-1]), atol=7), np.isclose(S, 0, atol=10**(-2)))
        #nullspace = np.isclose(np.log(S), np.log(S[-1]), atol=7)
        nullspace_rank = np.sum(nullspace)
        if nullspace_rank == 1:
            r_tilde_opt = np.transpose(np.atleast_2d(U[:, -1]))
            if np.isclose(r_tilde_opt[-1, 0], 0):
                test=1
            r_tilde_opt = r_tilde_opt / r_tilde_opt[-1, 0]
        elif nullspace_rank == 2:
            potential_vec = U[:, nullspace]
            non_zero = [not np.allclose(potential_vec[:9, i].flatten(), np.zeros((9))) for i in range(2)]
            if np.sum(non_zero)>1:
                return runtime, 0, 0, 0, 0, 0, False
            r_tilde_opt = potential_vec[:, non_zero]
            scalingfactor = 1. / (np.abs(np.linalg.det(np.reshape(r_tilde_opt[:9, :], (3, 3), 'F'))))**(1. / 3.) * np.sign(np.linalg.det(np.reshape(r_tilde_opt[:9, :], (3, 3), 'F')))
            r_tilde_opt = r_tilde_opt * scalingfactor
        else:
            return runtime, 0, 0, 0, 0, 0, False

        ext_rot = np.reshape(r_tilde_opt[:9, 0], (3,3), order='F')
        if np.linalg.norm(np.transpose(ext_rot) @ ext_rot - np.eye(3), 'fro') > 0.1:
            return runtime, 0, 0, 0, 0, 0, False
        ext_t = -np.matmul(np.matmul(np.linalg.inv(Q[:4, :4]), Q[:4, 4:]), r_tilde_opt)
        x_opt = np.vstack((ext_t, r_tilde_opt))
        estimated_extcal = np.vstack((np.hstack((ext_rot, ext_t[:3, :])), [0, 0, 0, 1]))

        #Check the duality gap
        primal_value = (np.transpose(x_opt).dot(Q)).dot(x_opt)
        duality_gap = primal_value - dual_value

        if (np.abs(duality_gap)/primal_value).item() > 0.001:
            return runtime, 0, 0, 0, 0, 0, False

        #Check the ground truth
        gt_value = (np.transpose(self.gt).dot(Q)).dot(self.gt)

        if verbose:
            print('Valid Solution Found')
        return runtime, gt_value, primal_value, duality_gap, x_opt, estimated_extcal, True

    def relax_solve(self, cons='RCH', verbose=True):
        if verbose:
            print("SDP Relaxation")
        #Construct the weight matrix
        Q = self.construct_Q()

        #Check Rank of Q_tt
        if np.linalg.matrix_rank(Q[:3, :3]) < 3:
            print('Q_{t,t} rank < 3 for data given! Exiting Solver.')
            return None

        Q_tilde = Q[4:, 4:] - np.matmul(np.matmul(Q[4:, :4], np.linalg.inv(Q[:4, :4])), Q[:4, 4:])


        # # Create a list of constraints
        A, b, p = cg.construct_A_and_b(constraints=cons)

        # # Define and solve the CVXPY problem.
        # # Create a symmetric matrix variable.
        X = cp.Variable((10,10), symmetric=True)

        # # Define Constraints
        constraints = [X >> 0]
        constraints += [
            cp.trace(A[i] @ X) == b[i] for i in range(p)
        ]

        # # Solve the Relaxed SDP
        prob = cp.Problem(cp.Minimize(cp.trace(Q_tilde @ X)),
                        constraints)
        #Initialize timer
        start_time = time.time()
        prob.solve(verbose=verbose, solver=cp.CVXOPT, abstol=10**(-7), refinement=1, kktsolver='robust')
        runtime = time.time() - start_time

        #Exit if solution was not found
        if prob.status in ["infeasible", "unbounded"]:
            return runtime, 0, 0, 0, 0, 0, False

        # Extract r and alpha from X
        X_value = X.value

        U, S, _ = np.linalg.svd(X_value)
        eigvec = np.logical_not(np.isclose(S, np.zeros(S.shape), atol=10**(-3)))
        rank = np.sum(eigvec)
        if rank ==1:
            r_tilde_opt = np.transpose(np.atleast_2d(U[:, 0] / U[-1, 0]))
        elif rank ==2:
            potential_solutions = U[:, eigvec]
            non_zero = [not np.allclose(potential_solutions[:9, i].flatten(), np.zeros((9))) for i in range(2)]
            if np.sum(non_zero)>1:
                return runtime, 0, 0, 0, 0, 0, False
            r_tilde_opt = potential_solutions[:, non_zero]
            scalingfactor = 1. / (np.abs(np.linalg.det(np.reshape(r_tilde_opt[:9, :], (3, 3), 'F')))) ** (
                        1. / 3.) * np.sign(np.linalg.det(np.reshape(r_tilde_opt[:9, :], (3, 3), 'F')))
            r_tilde_opt = r_tilde_opt * scalingfactor
        else:
            return runtime, 0, 0, 0, 0, 0, False

        ext_rot = np.reshape(r_tilde_opt[:9, 0], (3,3), order='F')
        if np.linalg.norm(np.transpose(ext_rot) @ ext_rot - np.eye(3), 'fro') > 10**(-3):
            return runtime, 0, 0, 0, 0, 0, False

        ext_t = -np.matmul(np.matmul(np.linalg.inv(Q[:4, :4]), Q[:4, 4:]), r_tilde_opt)
        x_opt = np.vstack((ext_t, r_tilde_opt))
        estimated_extcal = np.vstack((np.hstack((ext_rot, ext_t[:3, :])), [0, 0, 0, 1]))
        primal = np.transpose(x_opt) @ Q @ x_opt
        rel_gap = primal- prob.value

        if (np.abs(rel_gap)/primal).item() > 0.0001:
            return runtime, 0, 0, 0, 0, 0, False

        #Check the ground truth
        gt_value = (np.transpose(self.gt).dot(Q)).dot(self.gt)

        if verbose:
            print('Valid Solution Found')
        return runtime, gt_value, primal, rel_gap, x_opt, estimated_extcal, True

def main():
    #Name the relevant filenames
    data_filename = "unscaled_pose_data.mat"

    #Load the data
    T_vki, T_ci, extcal =load_data_1(data_filename)

    #Initialize solver
    my_solver = solver()

    #Convert inertial poses to relative
    T_v_rel = inertial_to_relative(T_vki)
    T_c_rel = inertial_to_relative(T_ci)

    #Choose a subset of the data
    num_poses_total = len(T_v_rel)
    num_poses = num_poses_total
    subset_indices = np.random.choice(num_poses_total, num_poses, replace=False).tolist()
    T_v_rel_sub = [T_v_rel[index] for index in subset_indices]
    T_c_rel_sub = [T_c_rel[index] for index in subset_indices]
    my_solver.set_T_v_rel(T_v_rel_sub)
    my_solver.set_T_c_rel(T_c_rel_sub, 2)
    my_solver.set_extcal(extcal)

    #Run Solver
    dual_time, dual_gt, dual_primal, dual_gap, dual_opt, dual_solution, opt_flag = my_solver.dual_solve(cons='R')
    rel_time, rel_gt, rel_primal, rel_gap, relax_opt, relax_solution, opt_flag = my_solver.relax_solve(cons='R')

    return

if __name__ == "__main__":
    main()
