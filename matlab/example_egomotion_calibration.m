%% Example usage of egomotion_calibration.m
close all; clear; clc;

% Create data matrices
n_points = 10;
% Noise levels
sigma_t = 0.01;
sigma = 0.01;
R1 = zeros(3,3,n_points);
R2 = zeros(3,3,n_points);
t1 = zeros(3,n_points);
t2 = zeros(3,n_points);
R = random_rot_matrix(3);
t = rand(3,1);
K = [R t; 0 0 0 1];
% Fill with random poses
for jdx=1:n_points
    R1(:,:,jdx) = random_rot_matrix(3);
    t1(:,jdx) = rand(3,1);
    Tjdx = [R1(:,:,jdx) t1(:,jdx); 0 0 0 1];
    T2jdx = K*Tjdx*inv(K);
    
    R2(:,:,jdx) = get_rotation_noise(sigma)*T2jdx(1:3,1:3);
    t2(:,jdx) = T2jdx(1:3,4) + randn(3,1)*sigma_t;
    R1(:,:,jdx) = get_rotation_noise(sigma)*R1(:,:,jdx);
    t1(:,jdx) = t1(:,jdx) + randn(3,1)*sigma_t;
end

% Solve
disp('Estimate:');
[R_cal, t_cal] = egomotion_calibration(R1,t1,R2,t2)
disp('Ground Truth:');
R
t