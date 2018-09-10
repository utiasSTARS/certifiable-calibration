

%% 'Minimal' observability problem analysis
% Don't use random noise - use magnitude of translation error added, and
% look at all 4 (row orthog, +handed, +col orthog, +handed + col orthog)
% possible configurations of constraints 
clear; close all; clc;

[r,p,w]=bot_quat_to_roll_pitch_yaw(bot_angle_axis_to_quat(pi/2, [0 1 0]));
trueCalib_euler = [0.06 0.03 0.04 r p w ];
n_points = 2;
tol = 1e-2;
% For reproducibility
rand('seed', 1);

%% Setup "minimal" problem
R1 = zeros(3,3,n_points);
R2 = zeros(3,3,n_points);
t1 = zeros(3,n_points);
t2 = zeros(3,n_points);

axis1 = [1; 0; 0];
alpha1 = pi/2;
axis2 = [0; 1; 0];
alpha2 = pi/2;
R1(:,:,1) = axis_angle(axis1, alpha1);
R1(:,:,2) = axis_angle(axis2, alpha2);
t1(:,1) = [1; 0; 0];
t1(:,2) = [0; 1; 0];

R = random_rot_matrix(3);
t = rand(3,1);
K = [R t; 0 0 0 1];

diff_cam1_dq = zeros(n_points,8);
diff_cam2_dq = zeros(n_points,8);

for jdx=1:n_points
    Tjdx = [R1(:,:,jdx) t1(:,jdx); 0 0 0 1];
    T2jdx = K*Tjdx*inv(K);
    R2(:,:,jdx) = T2jdx(1:3,1:3);
    t2(:,jdx) = T2jdx(1:3,4);
end

%% Run experiments
n_magnitudes = 4;
rot_magnitudes = linspace(pi/(n_magnitudes), pi, n_magnitudes);
n_res = 4;
u_vals = linspace(-1,1, n_res);
theta_vals = linspace(0, 2*pi, n_res+1);
theta_vals = theta_vals(1:end-1);
trans_magnitudes = [1 10 100 1000];

default_gap = zeros(n_res^2, n_res^2, n_magnitudes, n_magnitudes);
orthog_gap = zeros(n_res^2, n_res^2, n_magnitudes);
handed_gap = zeros(n_res^2, n_res^2, n_magnitudes);
orthog_handed_gap = zeros(n_res^2, n_res^2, n_magnitudes, n_magnitudes);

default_primal = zeros(n_res^2, n_res^2, n_magnitudes, n_magnitudes);
orthog_primal = zeros(n_res^2, n_res^2, n_magnitudes, n_magnitudes);
handed_primal = zeros(n_res^2, n_res^2, n_magnitudes, n_magnitudes);
orthog_handed_primal = zeros(n_res^2, n_res^2, n_magnitudes, n_magnitudes);

R2_in = R2;
t1_in = t1;
for idx=1:n_res
u = u_vals(idx);
for jdx=1:n_res
kdx = (idx-1)*n_res + jdx
theta = theta_vals(jdx);
axis_error = zeros(3,1);
axis_error(1) = sqrt(1-u^2)*cos(theta);
axis_error(2) = sqrt(1-u^2)*sin(theta);
axis_error(3) = u;
for rot_idx=1:n_magnitudes
    rot_mag = rot_magnitudes(rot_idx);
    R_err = axis_angle(axis_error, rot_mag);            
    R2_in(:,:,1) = R_err*R2(:,:,1);

    for idx_t=1:n_res
    u_t = u_vals(idx_t);
    for jdx_t=1:n_res
    kdx_t = (idx_t-1)*n_res + jdx_t
    theta_t = theta_vals(jdx_t);
    trans_error = zeros(3,1);
    trans_error(1) = sqrt(1-u^2)*cos(theta_t);
    trans_error(2) = sqrt(1-u^2)*sin(theta_t);
    trans_error(3) = u_t;
    
    for trans_idx=1:n_magnitudes
        trans_mag = trans_magnitudes(trans_idx);
        t_err = trans_error*trans_mag;
        t1_in(:,1) = t_err + t1(:,1);
        results = run_duality_experiment(R1,t1_in,R2_in,t2);

        default_gap(kdx, kdx_t, rot_idx, trans_idx) = results.default.gap;
        orthog_gap(kdx, kdx_t, rot_idx, trans_idx) = results.orthog.gap;
        handed_gap(kdx, kdx_t, rot_idx, trans_idx) = results.handed.gap;
        orthog_handed_gap(kdx, kdx_t, rot_idx, trans_idx) = results.orthog_handed.gap;

        default_primal(kdx, kdx_t, rot_idx, trans_idx) = results.default.primal_val;
        orthog_primal(kdx, kdx_t, rot_idx, trans_idx) = results.orthog.primal_val;
        handed_primal(kdx, kdx_t, rot_idx, trans_idx) = results.handed.primal_val;
        orthog_handed_primal(kdx, kdx_t, rot_idx, trans_idx) = results.orthog_handed.primal_val;
    end
    
end
end
end
end
end

% 
% a_err = [0; 0; 1];
% err_mag = -0.6*pi; %3*pi/4;
% R2(:,:,1) = R2(:,:,1)*axis_angle(a_err, err_mag);
% results = run_duality_experiment(R1,t1,R2,t2);

%% Post-process

default_count = sum(sum(abs(default_gap./default_primal) < tol, 1), 2);
orthog_count = sum(sum(abs(orthog_gap./orthog_primal) < tol, 1), 2);
handed_count = sum(sum(abs(handed_gap./handed_primal) < tol, 1), 2);
orthog_handed_count = sum(sum(abs(orthog_handed_gap./orthog_handed_primal) < tol, 1), 2);


save('../data/duality_gap/constraint_comparison_translation.mat');
