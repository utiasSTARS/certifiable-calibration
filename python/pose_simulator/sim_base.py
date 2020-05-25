import grd_gen as ggen
import traj_generator as tgen
import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from liegroups.numpy.se3 import SE3
from liegroups.numpy.so3 import SO3



class sim_base():
    def __init__(self, sim_info):
        #Store information from sim_info
        self.num_poses = sim_info['num_poses']
        self.dt = sim_info['dt']
        self.radius = sim_info['bounding_radius']
        self.traj_name = sim_info['traj_name']

        #Generate the ground and the path
        self.gnd_generator = self.gen_ground()
        self.tgen = self.gen_trajectory()
        self.extcal = self.gen_extrinsic()

        #Store the useful information
        self.t, self.T_vki_list, v, self.rot_vel = self.tgen.get_data() #time, Pose information, velocity in global reference frame, rotational velocity in sensor reference frame
        self.vel = np.empty(v.shape) #Velocity in pose sensor frame
        self.c_vel = np.empty(v.shape) #Velocity in camera sensor frame
        rcvkvk = self.extcal.inv().trans
        for i, pose in enumerate(self.T_vki_list):
            self.vel[:, i] = pose.rot.dot(v[:, i])
            self.c_vel[:, i] = (self.extcal.rot).dot(self.vel[:, i] + np.cross(self.rot_vel[:, i], rcvkvk)) #Camera velocity in camera frame

        #Create the camera data
        self.T_ci_list=[self.extcal.dot(pose) for pose in self.T_vki_list] #Camera Pose list
        self.c_rot_vel = np.transpose((self.extcal.rot).dot(np.transpose(self.rot_vel))) # Camera Rotational Velocity in the Camera reference frame

        return

    def gen_ground(self):
        gnd_generator = ggen.gnd_gen(None)
        return gnd_generator

    def gen_trajectory(self):
        traj_gen = tgen.traj_gen(self.num_poses, self.dt, self.radius, self.gnd_generator, None)
        return traj_gen

    def gen_extrinsic(self):
        #Generate the random translation vector
        rcvkvk = np.random.uniform(size=3)
        rot = np.random.uniform(-np.pi, np.pi, 3)
        Ccvk = SO3.exp(rot)
        ext_cal = SE3(rot=Ccvk, trans=-Ccvk.dot(rcvkvk))
        return ext_cal

    def SE3_to_npy(self, list):
        list_npy = [pose.as_matrix() for pose in list]
        array = np.array(list_npy)
        return array

    def sim_data_saver(self, data_list, filename):
        data = {}
        for (key, value) in data_list:
            if isinstance(value, list):
                if all(isinstance(item, SE3) for item in value):
                    value_npy=self.SE3_to_npy(value)
                    data[key] = value_npy
            elif isinstance(value, np.ndarray):
                data[key] = value
            elif isinstance(value, SE3):
                data[key] = value.as_matrix()
            else:
                print('Not SE(3) or ndarray trying to be saved!')
                break

        sio.savemat(filename, data)
        return

def visualize_world(sim):
    # Output figure of the world
    file_name = "./{}_sim_rand_data.pdf".format(sim.traj_name)
    plt.rc('text', usetex=True)
    fig = plt.figure()
    ax = Axes3D(fig)
    r_vi_i = np.empty((len(sim.T_vki_list), 3))
    r_ci_i = np.empty((len(sim.T_ci_list), 3))
    for i in range(len(sim.T_vki_list)):
        r_vi_i[i] = sim.T_vki_list[i].inv().trans
        r_ci_i[i] = sim.T_ci_list[i].inv().trans

    ax.plot3D(r_vi_i[:, 0], r_vi_i[:, 1], r_vi_i[:, 2], '-', \
                color='#609BCB', label=r'Egomotion Sensor Trajectory', linewidth=1.0)
    ax.plot3D(r_ci_i[:, 0], r_ci_i[:, 1], r_ci_i[:, 2], '-', \
                color='#E77554', label=r'Monocular Camera Trajectory', linewidth=1.0)
    ax.legend()
    ax.grid(True, which='both')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('z (m)')
    ax.set_alpha(None)
    # plt.title(r'Pose Sensor and Camera Trajectory')
    fig.savefig(file_name, bbox_inches='tight')
    return

def main():
    sim_setup = {
        "num_poses": 1000,  # Both num_poses and traj_names can be lists that correspond to different datasets
        "traj_name": 'pose_trajectory',
        "dt": 1. / 14,
        'bounding_radius': 50,
        "data_output_folder": ".",
    }
    sim = sim_base(sim_setup)
    visualize_world(sim)
    data = [('T_vki_list', sim.T_vki_list), ('T_ci_list', sim.T_ci_list), ('extcal', sim.extcal)]
    filename = 'unscaled_pose_data.mat'
    sim.sim_data_saver(data, filename)
    return

if __name__ == "__main__":
    main()