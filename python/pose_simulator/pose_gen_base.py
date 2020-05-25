#Base Class, Pose Class, and Odometry Class
import numpy as np
from liegroups.numpy.se3 import SE3
from liegroups.numpy.so3 import SO3

class PoseGenBase:
    def __init__(self, num_poses, semisphere_radius):
        self.num_poses = num_poses
        self.semisphere_radius = semisphere_radius
        self.T_vi_list = []

    def get_data(self):
        return self.T_vi_list

class PoseGen(PoseGenBase, object):
    def __init__(self, num_poses, semisphere_radius, angle_limits):
        super(PoseGen, self).__init__(num_poses, semisphere_radius)
        self.angle_limits = angle_limits
        self.gen_poses()
        return

    def gen_6dof_poses(self):
        # Sample random positions within the interior (defined by 50% radius distance) of the semi-sphere
        th = (np.pi / 2.) * np.random.rand(self.num_poses)
        phi = (2. * np.pi) * np.random.rand(self.num_poses)
        radii = 0.25 * self.semisphere_radius * np.random.rand(self.num_poses)

        x = radii * np.sin(th) * np.cos(phi)
        y = radii * np.sin(th) * np.sin(phi)
        z = radii * np.cos(th)

        # Select random angles for the rotation matrix such that is still faces roughly towards the semisphere
        (z_lim, y_lim, x_lim) = self.angle_limits
        angle_z = z_lim * np.random.rand(self.num_poses)
        angle_y = (y_lim * 2) * np.random.rand(self.num_poses) - y_lim
        angle_x = (x_lim * 2) * np.random.rand(self.num_poses) - x_lim
        return x, y, z, angle_x, angle_y, angle_z

    def gen_poses(self):
        x, y, z, angle_x, angle_y, angle_z = self.gen_6dof_poses()
        for i in range(self.num_poses):
            r_vi_i = np.array([x[i], y[i], z[i]])
            C_vi = SO3.rotx(angle_x[i]).dot(SO3.roty(angle_y[i])).dot(SO3.rotz(angle_z[i]))
            self.T_vi_list.append(SE3(rot=C_vi, trans=-1 * (C_vi.dot(r_vi_i))))
        return

class OdomGen(PoseGen, object):
    def __init__(self, num_poses, semisphere_radius, angle_limits):
        super(OdomGen, self).__init__(num_poses, semisphere_radius, angle_limits)
        self.gen_odom()
        return

    def gen_odom(self):
        temp_list = []
        for pose in self.T_vi_list:
            angle = 3. * (3.1415 / 180.)
            rand_vec = np.random.randn(3)
            rand_vec = rand_vec / np.linalg.norm(rand_vec)
            C_v2v1 = SO3.exp(angle * rand_vec)
            r_v2v1_v1 = 0.1 * np.random.randn(3)
            T_v2v1 = SE3(rot=C_v2v1, trans=-1 * (C_v2v1.dot(r_v2v1_v1)))

            T_v2i = T_v2v1.dot(pose)

            temp_list.append(pose)
            temp_list.append(T_v2i)
        self.T_vi_list = temp_list
        return