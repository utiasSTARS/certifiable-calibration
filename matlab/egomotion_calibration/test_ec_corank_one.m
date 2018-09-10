%% test co-rank 1

clear; close all; clc;


%% Zero noise case
n_points = 20; %10;
tol = 1e-4;

R1 = zeros(3,3,n_points);
R2 = zeros(3,3,n_points);
t1 = zeros(3,n_points);
t2 = zeros(3,n_points);
R = random_rot_matrix(3);
t = rand(3,1);
K = [R t; 0 0 0 1];

diff_cam1_dq = zeros(n_points,8);
diff_cam2_dq = zeros(n_points,8);

 % assign an equally weighted covariance for all measurements in the lie algebra
cov1_dq = eye(6*size(diff_cam1_dq,1)) * 1e-5;
cov2_dq = eye(6*size(diff_cam1_dq,1)) * 1e-5;

axis = rand(3,1);
axis = axis/norm(axis);

for jdx=1:n_points

    axis = rand(3,1);
    axis = axis/norm(axis);
    R1(:,:,jdx) = axis_angle(axis, rand*2*pi);

    t1(:,jdx) = rand(3,1);
    diff_cam1_dq(jdx,:) = EulerStateToDualQuat([t1(:,jdx).',rotm2euler(R1(:,:,jdx))]);
    
    Tjdx = [R1(:,:,jdx) t1(:,jdx); 0 0 0 1];
    T2jdx = K*Tjdx*inv(K);
    R2(:,:,jdx) = T2jdx(1:3,1:3);
    t2(:,jdx) = T2jdx(1:3,4);
    diff_cam2_dq(jdx,:) = EulerStateToDualQuat([t2(:,jdx).',rotm2euler(R2(:,:,jdx))]);
end

M = ec_get_data_matrix(R1,t1,R2,t2);
[A,b,c] = ec_get_affine_cost_matrices(R1,t1,R2,t2);
% Get Schur complement
Mrr = M(4:end, 4:end);
Mtt = M(1:3, 1:3);
Mrt = M(1:3, 4:end);
Q = Mrr - Mrt.'*inv(Mtt)*Mrt;

%% Step 1: Test Strict convexity of non-homogenized cost function
eigs_affine = eig(A)

%% Step 2: Test Corank-one of Lagrangian Hessian

%% Compare all strengthening options here
% TODO: Fix interface to independently consider boolean options for
% strengthening (so that we can test handedness without orthog strength)
% Also do ACQ??? Or is the "variety radical"?
[dual_val, dual_sol] = ec_dual_solver(Q, true, true);
eigs_strengthened = eig(dual_sol.H);
dual_val

R_dual = ec_extract_primal(dual_sol);
t_dual = ec_get_t_opt(R_dual, Mtt, Mrt);
