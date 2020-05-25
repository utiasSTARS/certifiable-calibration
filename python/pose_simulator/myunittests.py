import numpy as np
import sim_base as sb


sim_setup = {
    "num_poses": 1000,  #Both num_poses and traj_names can be lists that correspond to different datasets
    "traj_name": 'pose_trajectory',
    "dt" : 1./14,
    'bounding_radius': 50,
    "data_output_folder": ".",
}

def check_velocity_pair(vel_1, vel_2, rot_vel2, ext_cal):
    r_sc_c = ext_cal.trans
    C_sc = ext_cal.inv().rot
    error = np.empty(vel_1.shape[1])
    for i in range(len(vel_2)):
        error[i] = np.linalg.norm(vel_1[:, i] - C_sc.dot(vel_2[:, i] + np.cross(rot_vel2[:, i], r_sc_c)))
    return

def check_velocity_integration(dt, pose_list, vel): #Verify the velocity integration is close to the translation value
    #Unpack the pose information
    r_si_i = np.empty((3, len(pose_list)))
    v_c_i = np.empty((3, len(pose_list)))
    for i, pose in enumerate(pose_list):
        r_si_i[:, i] = pose.inv().trans
        v_c_i[:, i] = np.transpose(pose.inv().rot.dot(vel[:, i]))
    d_c_i_true = np.diff(r_si_i)
    d_c_i = dt * v_c_i
    error = np.linalg.norm(d_c_i_true - d_c_i[:, :len(pose_list)-1], axis=0)
    return

def check_poses(pose_list_1, pose_list_2, ext_cal):
    rel_motion_1 = [ext_cal.dot(pose_list_1[i+1].dot(pose_list_1[i].inv())) for i in range(len(pose_list_1)-1)]
    rel_motion_2 = [(pose_list_2[i+1].dot(pose_list_2[i].inv())).dot(ext_cal) for i in range(len(pose_list_1)-1)]
    error = [rel_motion_1[i].dot(rel_motion_2[i].inv()) for i in range(len(rel_motion_1))]
    return

def main():
    sim = sb.sim_base(sim_setup)
    check_velocity_integration(sim_setup["dt"], sim.T_ci_list, sim.c_vel)
    check_poses(sim.T_vki_list, sim.T_ci_list, sim.extcal)
    check_velocity_pair(sim.vel, sim.c_vel, sim.c_rot_vel, sim.extcal)
    return

if __name__ == "__main__":
    main()