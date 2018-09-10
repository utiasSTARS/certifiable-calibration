%% Test ec_dual_solver
clear; close all; clc;


%% Zero noise case
n_points = 10;
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

for jdx=1:n_points
    R1(:,:,jdx) = random_rot_matrix(3);
    t1(:,jdx) = rand(3,1);
    diff_cam1_dq(jdx,:) = EulerStateToDualQuat([t1(:,jdx).',rotm2euler(R1(:,:,jdx))]);
    
    Tjdx = [R1(:,:,jdx) t1(:,jdx); 0 0 0 1];
    T2jdx = K*Tjdx*inv(K);
    R2(:,:,jdx) = T2jdx(1:3,1:3);
    t2(:,jdx) = T2jdx(1:3,4);
    diff_cam2_dq(jdx,:) = EulerStateToDualQuat([t2(:,jdx).',rotm2euler(R2(:,:,jdx))]);
end

[estimatedCalibration_dq, optimOutput] = Calibrate3d_DualQuat3(diff_cam2_dq, cov2_dq, diff_cam1_dq, cov1_dq, [], 'zeros', 0);
estimatedCalibration_euler = DualQuatToEulerState(estimatedCalibration_dq);

M = ec_get_data_matrix(R1,t1,R2,t2);
% Get Schur complement
Mrr = M(4:end, 4:end);
Mtt = M(1:3, 1:3);
Mrt = M(1:3, 4:end);
Q = Mrr - Mrt.'*inv(Mtt)*Mrt;
[dual_val, dual_sol] = ec_dual_solver(Q, true, true);
R_dual = ec_extract_primal(dual_sol);
t_dual = ec_get_t_opt(R_dual, Mtt, Mrt);

disp(estimatedCalibration_euler);
disp([t_dual.' , rotm2euler(R_dual)]);
disp([t.' , rotm2euler(R)]);

assert(all(all(abs(R-R_dual) < tol)), 'Rotation not close enough!');
assert(all(all(abs(t-t_dual) < tol)), 'Translation not close enough!');

%% Noisy case

n_points = 100;
sigmas = linspace(0.0, 1.0, 10);
sigmas_t = sigmas;
gap_default = zeros(size(sigmas));
gap_orthog = zeros(size(sigmas));
gap_handed = zeros(size(sigmas));


for idx=1:length(sigmas)
    sigma = sigmas(idx);
    sigma_t = sigmas_t(idx);
    R1 = zeros(3,3,n_points);
    R2 = zeros(3,3,n_points);
    t1 = zeros(3,n_points);
    t2 = zeros(3,n_points);
    R = random_rot_matrix(3);
    t = rand(3,1);
    K = [R t; 0 0 0 1];

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

    M = ec_get_data_matrix(R1,t1,R2,t2);
    % Get Schur complement
    Mrr = M(4:end, 4:end);
    Mtt = M(1:3, 1:3);
    Mrt = M(1:3, 4:end);
    Q = Mrr - Mrt.'*inv(Mtt)*Mrt;
    
%     [dual_val, dual_sol] = ec_dual_solver(Q, false, false);
%     R_dual = ec_extract_primal(dual_sol);
%     t_dual = ec_get_t_opt(R_dual, Mtt, Mrt);
%     x_dual = [vec(t_dual); vec(R_dual); 1];
%     primal_val = x_dual.'*M*x_dual;
%     gap_default(idx) = primal_val - dual_val;
    
%     [dual_val, dual_sol] = ec_dual_solver(Q, true, false);
%     R_dual = ec_extract_primal(dual_sol);
%     t_dual = ec_get_t_opt(R_dual, Mtt, Mrt);
%     x_dual = [vec(t_dual); vec(R_dual); 1];
%     primal_val = x_dual.'*M*x_dual;
%     gap_orthog(idx) = primal_val - dual_val;
    
    [dual_val, dual_sol] = ec_dual_solver(Q, true, true);
    R_dual = ec_extract_primal(dual_sol);
    t_dual = ec_get_t_opt(R_dual, Mtt, Mrt);
    x_dual = [vec(t_dual); vec(R_dual); 1];
    primal_val = x_dual.'*M*x_dual;
    gap_handed(idx) = primal_val - dual_val;
    
%     assert(abs(dual_val-primal_val) < tol, 'Strong duality not holding');
    
    
end


%% Visualize gap
figure; hold on; grid on;
plot(sigmas, gap_default, 'r-o');
plot(sigmas, gap_orthog, 'g--s');
plot(sigmas, gap_handed, 'b-.d');
legend('Default', 'Orthogonality Str.', 'Handed Str.');
xlabel('Sigma');
ylabel('Duality Gap (p-d)');


