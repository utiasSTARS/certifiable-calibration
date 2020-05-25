#Trajectory Class
from pose_gen_base import PoseGenBase
import numpy as np
from liegroups.numpy.se3 import SE3
from liegroups.numpy.so3 import SO3

def normalize(x):
    return x / np.linalg.norm(x)

class traj_gen(PoseGenBase, object):
    def __init__(self, num_poses, dt, radius,  gnd, constants):
        super(traj_gen, self).__init__(num_poses, radius)
        #Calculate the simulation time
        self.t = np.arange(0, self.num_poses) * dt

        #Variables to store info
        self.rot_vel = np.empty((3, num_poses))
        self.vel = np.empty((3, num_poses))

        #Parameters to define the surface
        if constants is None:
            self.alpha = [0.5 * self.semisphere_radius, 0.1 * self.semisphere_radius, np.random.rand(),
                          (0.5 + 0.1 * np.random.rand()),
                          (0.5 + 0.1 * np.random.rand())]
        else:
            self.alpha = constants
        self.beta = gnd.constants

        #Construct the data
        self.calc_motion_data(gnd)
        self.print_stats()
        return

    def calc_r(self):
        r = self.alpha[0] + self.alpha[1] * np.sin(self.alpha[2] * self.t)
        return r

    def calc_dr_dt(self):
        dr_dt = self.alpha[1] * self.alpha[2] * np.cos(self.alpha[2] * self.t)
        return dr_dt

    def calc_ddr_dtdt(self):
        ddr_dtdt = - self.alpha[1] * self.alpha[2] ** 2 * np.sin(self.alpha[2] * self.t)
        return ddr_dtdt

    def calc_x(self, r):
        x = r * np.cos(self.alpha[3] * self.t)
        return x

    def calc_dx_dt(self, r, dr_dt):
        dx_dt = dr_dt * np.cos(self.alpha[3] * self.t) - self.alpha[3] * r * np.sin(self.alpha[3] * self.t)
        return dx_dt

    def calc_ddx_dtdt(self, r, dr_dt, ddr_dtdt):
        ddx_dtdt = ddr_dtdt * np.cos(self.alpha[3] * self.t) - 2 * self.alpha[3] * dr_dt * np.sin(self.alpha[3] * self.t) - self.alpha[3] ** 2 * r * np.cos(self.alpha[3] * self.t)
        return ddx_dtdt

    def calc_y(self, r):
        y = r * np.sin(self.alpha[4] * self.t)
        return y

    def calc_dy_dt(self, r, dr_dt):
        dy_dt = dr_dt * np.sin(self.alpha[4] * self.t) + r * self.alpha[4] * np.cos(self.alpha[4] * self.t)
        return dy_dt

    def calc_ddy_dtdt(self, r, dr_dt, ddr_dtdt):
        ddy_dtdt = ddr_dtdt * np.sin(self.alpha[4] * self.t) + 2 * self.alpha[4] * dr_dt * np.cos(self.alpha[4] * self.t) - self.alpha[4] ** 2 * r * np.sin(self.alpha[4] * self.t)
        return ddy_dtdt

    def calc_z(self, x, y, gnd):
        z = gnd.calc_z(x, y)
        return z

    def calc_dz_dx(self, x, y, gnd):
        dz_dx = gnd.calc_dz_dx(x, y)
        return dz_dx

    def calc_ddz_dxdx(self, x, y, gnd):
        ddz_dxdx = gnd.calc_ddz_dxdx(x, y)
        return ddz_dxdx

    def calc_dz_dy(self, x, y, gnd):
        dz_dy = gnd.calc_dz_dy(x, y)
        return dz_dy

    def calc_ddz_dydy(self, x, y, gnd):
        ddz_dydy = gnd.calc_ddz_dydy(x, y)
        return ddz_dydy

    def calc_ddz_dxdy(self, x, y, gnd):
        ddz_dxdy = gnd.calc_ddz_dxdy(x, y)
        return ddz_dxdy

    def calc_dz_dt(self, dx_dt, dy_dt, dz_dx, dz_dy):
        dz_dt = dz_dx * dx_dt + dz_dy * dy_dt
        return dz_dt

    def calc_ddz_dxdt(self, dx_dt, dy_dt, ddz_dxdx, ddz_dxdy):
        ddz_dxdt = ddz_dxdx * dx_dt + ddz_dxdy * dy_dt
        return ddz_dxdt

    def calc_ddz_dydt(self, dx_dt, dy_dt, ddz_dydy, ddz_dxdy):
        ddz_dydt = ddz_dxdy * dx_dt + ddz_dydy * dy_dt
        return ddz_dydt

    def calc_ddz_dtdt(self, dz_dx, ddz_dxdx, ddz_dxdy, dz_dy, ddz_dydy, dx_dt, ddx_dtdt, dy_dt, ddy_dtdt):
        ddz_dtdt = ddz_dxdx * np.square(dx_dt) + 2 * ddz_dxdy * dy_dt * dx_dt + ddz_dydy * np.square(dy_dt) + dz_dx * ddx_dtdt + dz_dy * ddy_dtdt
        return ddz_dtdt

    def calc_comm_path_values(self, gnd):
        #Zeroth Derivatives
        r = self.calc_r()
        dr_dt = self.calc_dr_dt()
        ddr_dtdt = self.calc_ddr_dtdt()
        x = self.calc_x(r)
        y = self.calc_y(r)
        z = self.calc_z(x, y, gnd)
        return r, dr_dt, ddr_dtdt, x, y, z

    def calc_translation(self, x, y, z):
        trans = np.vstack((x, y, z))
        return trans

    def calc_vectors(self, r, dr_dt, ddr_dtdt, x, y, gnd):
        # Calculate the relavent derivatives
        # First Spatial Derivatives
        dz_dx = self.calc_dz_dx(x, y, gnd)
        dz_dy = self.calc_dz_dy(x, y, gnd)

        # First Time Derivatives
        dx_dt = self.calc_dx_dt(r, dr_dt)
        dy_dt = self.calc_dy_dt(r, dr_dt)
        dz_dt = self.calc_dz_dt(dx_dt, dy_dt, dz_dx, dz_dy)

        # Second Spatial Derivatives
        ddz_dxdx = self.calc_ddz_dxdx(x, y, gnd)
        ddz_dxdy = self.calc_ddz_dxdy(x, y, gnd)
        ddz_dydy = self.calc_ddz_dydy(x, y, gnd)

        #Calculate vectors
        tangents = self.tangent_vecs(dx_dt, dy_dt, dz_dt)
        normals = self.normal_vecs(dz_dx, dz_dy)
        atangents = self.tangent_acceleration_vecs(r, dr_dt, ddr_dtdt, dz_dx, ddz_dxdx, ddz_dxdy, dz_dy, ddz_dydy, dx_dt, dy_dt)
        anormals = self.normal_acceleration_vecs(dx_dt, dy_dt, ddz_dxdx, ddz_dxdy, ddz_dydy)
        return tangents, normals, atangents, anormals

    def normal_vecs(self, dz_dx, dz_dy):
        n = np.vstack((-dz_dx, -dz_dy, np.ones(self.t.shape[0])))
        return n

    def tangent_vecs(self, dx_dt, dy_dt, dz_dt):
        v = np.vstack((dx_dt, dy_dt, dz_dt))
        return v

    def tangent_acceleration_vecs(self, r, dr_dt, ddr_dtdt, dz_dx, ddz_dxdx, ddz_dxdy, dz_dy, ddz_dydy, dx_dt, dy_dt):
        ddx_dtdt = self.calc_ddx_dtdt(r, dr_dt, ddr_dtdt)
        ddy_dtdt = self.calc_ddy_dtdt(r, dr_dt, ddr_dtdt)
        ddz_dtdt = self.calc_ddz_dtdt(dz_dx, ddz_dxdx, ddz_dxdy, dz_dy, ddz_dydy, dx_dt, ddx_dtdt, dy_dt, ddy_dtdt)
        dv_dt = np.vstack((ddx_dtdt, ddy_dtdt, ddz_dtdt))
        return dv_dt

    def normal_acceleration_vecs(self, dx_dt, dy_dt, ddz_dxdx, ddz_dxdy, ddz_dydy):
        ddz_dxdt = self.calc_ddz_dxdt(dx_dt, dy_dt, ddz_dxdx, ddz_dxdy)
        ddz_dydt = self.calc_ddz_dydt(dx_dt, dy_dt, ddz_dydy, ddz_dxdy)
        dn_dt = np.vstack((-ddz_dxdt, -ddz_dydt, np.zeros(ddz_dxdt.shape)))
        return dn_dt

    def calc_SO3_matrix(self, e_1, e_2, e_3):
        SO3_matrix = np.stack((e_1, e_2, e_3), axis=1)
        return SO3_matrix

    def nmat_dot_nvec(self, mats, vecs):
        #Multiply N matrices with N associated vectors
        new_vecs = np.sum(mats * np.tile(vecs[None, :, :], (3, 1, 1)), axis=1)
        return new_vecs

    def unit_vector_der(self, v, dv_dt):
        # Calculate the time derivative of the tangent or normal unit vector
        v_mag = np.tile(np.linalg.norm(v, axis=0)[None, None, :], (3, 3, 1))
        id_tensor = np.tile(np.identity(3)[:, :, None], (1, 1, v.shape[1]))
        de_i_dv = id_tensor / v_mag - (v[:, None, :] * v[None, :, :]) / np.power(v_mag, 3)
        de_i_dt = self.nmat_dot_nvec(de_i_dv, dv_dt)
        return de_i_dt

    def binormal_unit_vector_der(self, e_1, e_2, de_1_dt, de_2_dt):
        #Calculate the derivative of the binormal unit vector
        de_3_dt = np.cross(de_1_dt, e_2, axisa=0, axisb=0) + np.cross(e_1, de_2_dt, axisa=0, axisb=0)
        return de_3_dt.T

    def calc_rot_vel(self, C_vi, e_1, e_2, e_3, de_1_dt, de_2_dt, de_3_dt):
        #Initialize Variables
        w_i = np.empty(e_1.shape)

        #Stack the basis vectors to solve for the rotation
        omega_skew = de_1_dt[None, :, :] * e_1[:, None, :] + de_2_dt[None, :, :] * e_2[:, None, :] + de_3_dt[None, :, :] * e_3[:, None, :]

        #Undo the skew symmetric operation and rotate into the vehicle reference frame
        w_i[0, :] = omega_skew[2, 1, :]
        w_i[1, :] = omega_skew[0, 2, :]
        w_i[2, :] = omega_skew[1, 0, :]
        w_v = self.nmat_dot_nvec(C_vi, w_i) #I don't know why this negative sign is needed? But it seems to produce the correct rot vel

        return w_v

    def calc_motion_data(self, gnd):
        #Calculate the required vectors
        r, dr_dt, ddr_dtdt, x, y, z = self.calc_comm_path_values(gnd)
        v, n, dv_dt, dn_dt = self.calc_vectors(r, dr_dt, ddr_dtdt, x, y, gnd)

        #Calculate the unit vectors
        e_1 = v / np.tile(np.linalg.norm(v, axis=0)[None, :], (3, 1))
        e_3 = n / np.tile(np.linalg.norm(n, axis=0)[None, :], (3, 1))
        e_2 = np.cross(e_3, e_1, axisa=0, axisb=0).T

        #Calculate the unit vector derivatives
        de_1_dt = self.unit_vector_der(v, dv_dt)
        de_3_dt = self.unit_vector_der(n, dn_dt)
        de_2_dt = self.binormal_unit_vector_der(e_3, e_1, de_3_dt, de_1_dt)

        #Calculate the translations
        trans = np.vstack((x, y, z))

        #Calculate the SO3 matrices
        SO3_matrices = self.calc_SO3_matrix(e_1, e_2, e_3)

        #Calculate the SE3 matrices
        for i in range(self.num_poses):
            C_iv = SO3.from_matrix(SO3_matrices[:, :, i])
            C_vi = C_iv.inv()
            r_ivk_vk = -C_vi.dot(trans[:, i])
            self.T_vi_list.append(SE3(rot=C_vi, trans=r_ivk_vk))

        #Store the velocity
        self.vel= v

        #Calculate the rotational velcities
        self.rot_vel = self.calc_rot_vel(np.transpose(SO3_matrices, (1, 0, 2)), e_1, e_2, e_3, de_1_dt, de_2_dt, de_3_dt)
        return

    def print_stats(self):
        #Print Speed statistics
        speeds = np.linalg.norm(self.vel, axis=0)
        avg_speed = speeds.mean()
        vel_var = np.mean(np.diff(self.vel, axis= 1), axis = 0)
        print('Max Speed: {:3.3f}'.format(np.max(speeds)))
        print('Average speed: {:3.3f}'.format(avg_speed))
        print('Velocity Mean Change: ({},{},{})'.format(vel_var[0], vel_var[1], vel_var[2]))

        #Print Rotational speed statistics
        rot_speed = np.linalg.norm(self.rot_vel, axis=0)
        rot_var_max = np.max(np.sum(np.diff(self.rot_vel, axis = 1), axis=0), axis=0)
        rot_var = np.mean(np.diff(self.rot_vel, axis = 1), axis=0)
        print('Max Rotational Speed: {:3.3f}'.format(np.max(rot_speed)))
        print('Average speed: {:3.3f}'.format(rot_speed.mean()))
        print('Rotational Velocity Mean Change: ({}, {}, {})'.format(rot_var[0], rot_var[1], rot_var[2]))
        print('Max Rotational Velocity Change: {}'.format(rot_var_max))
        return

    def get_data(self):
        return self.t, self.T_vi_list, self.vel, self.rot_vel
