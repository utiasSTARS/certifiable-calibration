import numpy as np
import math

class schur_solver:
    def __init__(self):
        T_c_rel = []
        R_v_rel = []
        return

    def set_vals(self, T_v_rel, T_c_rel):
        self.T_v_rel = T_c_rel #Reversed should be fixed before we open source the code
        self.T_c_rel = T_v_rel
        return True

    def rotation_projection(self, R):
        u, s, vh = np.linalg.svd(R, full_matrices=True)

        return u.dot(vh)

    def calculate_rotation(self):
        Q = self.construct_Q()

        u, s, vh = np.linalg.svd(Q, full_matrices=True)

        Rx = np.reshape(u[4:-1, -2], (3, 3))
        Rx = self.rotation_projection(Rx)
        Rx = Rx * np.sign(np.linalg.det(Rx))
        Rx = Rx.transpose()
        return Rx

    def solve(self):
        Rx = self.calculate_rotation()

        r = np.zeros((10, 1))
        r[:9] = np.reshape(Rx.transpose(), (9, 1))
        r[9] = 0

        Q = self.construct_Q()
        Qta = Q[:4, :4]
        Qtar = Q[:4, 4:]
        v = np.resize(-np.linalg.inv(Qta).dot(Qtar.dot(r)), (4, 1))

        tx = np.resize(v[:3], (3))
        scale = v[3]

        return [Rx, tx, scale]  #The inverse of scale is calculated, needs to be fixed

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

    def construct_Q(self):
        # Construct the rotation error cost
        M_R_list = self.construct_M_R()
        Q_R = self.construct_Q_i(M_R_list)

        #Construct the translation error cost
        M_t_list = self.construct_M_t()
        Q_t = self.construct_Q_i(M_t_list)

        Q = Q_t + Q_R
        Q = (np.transpose(Q) + Q)/2 #Make sure that Q is symmetric
        return Q
