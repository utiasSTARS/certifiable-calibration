import time
import numpy as np
import math
from utility.utils import load_data_1, inertial_to_relative
from solver.schur_solver import schur_solver

class solver():
    def __init__(self):
        #T_c is treated as A, T_v is treated as B
        self.T_c_rel = []
        self.T_v_rel = []
        self.extcal = []

        return

    def set_vals(self, T_v_rel, T_c_rel, extcal):
        self.T_c_rel = T_c_rel #Reversed should be fixed before we open source the code
        self.T_v_rel = T_v_rel
        self.extcal = extcal
        return True

    def rotation_projection(self, R):
        u, s, vh = np.linalg.svd(R, full_matrices=True)

        return u.dot(vh)

    def calculate_rotation_direct(self):
        count = np.shape(self.T_v_rel)[0]

        A = np.zeros((9 * count, 9))
        for i in range(count):
            Rai = self.T_c_rel[i].as_matrix()[:3, :3]
            Rbi = self.T_v_rel[i].as_matrix()[:3, :3]

            A[9*i:9*(i+1), :] = np.eye(9) - np.kron(Rai, Rbi)

        u, s, vh = np.linalg.svd(A, full_matrices=True)
        V = np.reshape(vh[-1,:].transpose(), (3, 3))
        Vdeterminant = np.linalg.det(V)

        R = V * np.sign(Vdeterminant) / math.pow(abs(Vdeterminant), 1/3)
        R = self.rotation_projection(R)

        return R

    def calculate_translation_scale_direct(self, Rx):
        count = np.shape(self.T_v_rel)[0]
        rx = np.resize(Rx, (9, 1))

        A = np.zeros((3*count, 4))
        b = np.zeros((3*count, 1))
        for i in range(count):
            Rai = self.T_c_rel[i].as_matrix()[:3, :3]
            Rbi = self.T_v_rel[i].as_matrix()[:3, :3]
            uai = self.T_c_rel[i].as_matrix()[:3, 3]
            tbi = self.T_v_rel[i].as_matrix()[:3, 3]

            A[3*i:3*(i+1), :3] = np.eye(3) - Rai
            A[3*i:3*(i+1), 3] = -uai
            b[3*i:3*(i+1)] = -np.kron(np.eye(3), tbi.transpose()).dot(rx)

        v = np.linalg.inv(A.transpose().dot(A)).dot(A.transpose()).dot(b)

        tx = v[:3].transpose()
        scale = v[3]
        return [tx, scale]

    def calculate_schur(self):
        s = schur_solver()
        s.set_vals(self.T_c_rel, self.T_v_rel)

        [Rx, tx, scale] = s.solve()

        return [Rx, tx, scale]

    def solve(self, schur=True):
        start_time = time.time()

        if schur:
            Rx, tx, scale = self.calculate_schur()
        else:
            Rx = self.calculate_rotation_direct()
            tx, scale = self.calculate_translation_scale_direct(Rx)

        runtime = time.time() - start_time
        X = np.eye(4)
        X[:3, :3] = Rx
        X[:3, 3] = tx

        return X, scale, runtime

    def report(self, schur=True):
        X, scale, runtime = self.solve(schur)

        print("Ground truth:")
        print(self.extcal)
        print("Solution:")
        print(X)
        print("Calculated scale:")
        print(scale)

        return True

def main():
    #Address of the input file
    data_filename = "unscaled_pose_data.mat"

    #Load the data
    T_v, T_c, extcal = load_data_1(data_filename)

    #Convert inertial poses to relative
    T_v_rel = inertial_to_relative(T_v)
    T_c_rel = inertial_to_relative(T_c)

    #Choose a subset of the data
    num_poses_total = len(T_v_rel)
    num_poses = 20

    subset_indices = np.random.choice(num_poses_total, num_poses, replace=False).tolist()
    T_v_rel_sub = [T_v_rel[index] for index in subset_indices]
    T_c_rel_sub = [T_c_rel[index] for index in subset_indices]

    #Initialize solver
    my_solver = solver()
    my_solver.set_vals(T_v_rel_sub, T_c_rel_sub, extcal)

    my_solver.report()

    return True


if __name__ == "__main__":
    main()
