import scipy.io as sio
import numpy as np
from liegroups import SE3, SO3

def load_data_1(data_filename):
    #Data loader for Emmett's Synethetic datasets
    data_dict = sio.loadmat(data_filename)
    pose_array = data_dict['T_vki_list']
    T_vki = [SE3.from_matrix(pose_array[i, :, :]) for i in range(pose_array.shape[0])]
    pose_array = data_dict['T_ci_list']
    T_ci = [SE3.from_matrix(pose_array[i, :, :]) for i in range(pose_array.shape[0])]
    extcal = data_dict['extcal']
    return T_vki, T_ci, extcal

def load_brookshire_data(data_filename):
    #Data loader for Emmett's Synethetic datasets
    data_dict = sio.loadmat(data_filename)
    pose_array = data_dict['T_v_rel_list']
    T_v_rel = [SE3.from_matrix(pose_array[:, :, i]) for i in range(pose_array.shape[2])]
    pose_array = data_dict['T_c_rel_list']
    T_c_rel = [SE3.from_matrix(pose_array[:, :, i]) for i in range(pose_array.shape[2])]
    extcal = data_dict['extcal']
    return T_v_rel, T_c_rel, extcal

def inertial_to_relative(pose_list, offset=1):
    #Changes SE3 poses from relative to intertial to relative vs timestamps
    rel_pose_list = [pose_list[i].dot(pose_list[i + offset].inv()) for i in range(0, len(pose_list) - offset, offset)] #This makes the poses T_{s_is_i+1}
    return rel_pose_list

def add_noise(pose_list, trans_noise, rot_noise):
    ''' Add Noise to the pose measurements as a percentage of the translation and rotation magnitude'''

    #Create the translation noise vectors
    trans_noise_sigma = trans_noise * np.linalg.norm(np.array([pose.trans  for pose in pose_list]), axis=1)
    trans_noise_vec = [np.random.randn((3))*trans_noise_sigma[i] for i in range(len(pose_list))]

    #Create the rotation noise vectors
    rot_noise_sigma = rot_noise * np.linalg.norm(np.array([pose.rot.log()  for pose in pose_list]), axis=1)
    rot_noise_unit_vecs = np.random.rand(3, len(pose_list))
    rot_noise_unit_vecs=  rot_noise_unit_vecs/ np.linalg.norm(rot_noise_unit_vecs, axis=0)
    rot_noise_vec = [rot_noise_sigma[i] * np.random.rand() * rot_noise_unit_vecs[:, i]  for i in range(len(pose_list))]

    #Add the noise to the poses
    noisy_poses = [SE3(rot=pose.rot.dot(SO3.exp(rot_noise_vec[i])), trans=pose.trans + trans_noise_vec[i]) for i, pose in enumerate(pose_list)]
    return noisy_poses
