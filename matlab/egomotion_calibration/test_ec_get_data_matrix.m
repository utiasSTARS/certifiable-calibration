%% test ec_get_data_matrix
clear; close all; clc;

%% Test Frobenius formulation
R1 = zeros(3,3,1);
R2 = zeros(size(R1));
R1(:,:,1) = eye(3);
R2(:,:,1) = eye(3);
t1 = zeros(3,1);
t2 = zeros(3,1);
c = primal_frobenius(eye(3), [0;0;0], R1, R2, t1, t2);
assert(c == 0);

%% Test cost function from quadratic form and affine form

n_trials = 100;
n_points = 100;
n_inner_trials = 10;
tol = 1e-10;
for idx=1:n_trials
    R1 = zeros(3,3,n_points);
    R2 = zeros(3,3,n_points);
    t1 = zeros(3,n_points);
    t2 = zeros(3,n_points);
    for jdx=1:n_points
        R1(:,:,jdx) = random_rot_matrix(3);
        R2(:,:,jdx) = random_rot_matrix(3);
        t1(:,jdx) = rand(3,1);
        t2(:,jdx) = rand(3,1);
    end
    M = ec_get_data_matrix(R1,t1,R2,t2);
    [A,b,c] = ec_get_affine_cost_matrices(R1,t1,R2,t2);
    for jdx=1:n_inner_trials
        R = random_rot_matrix(3);
        t = rand(3,1);
        x = [t; vec(R); 1];
        cf = primal_frobenius(R,t,R1,R2,t1,t2);
        cq = x.'*M*x;
        ca = x(1:12).'*A*x(1:12) + b.'*x(1:12) + c;
%         cf 
%         cq
%         cf-cq
        assert(abs(cf - cq) < tol, 'Homog. cost function incorrect!');
        assert(abs(cf - ca) < tol, 'Affine cost function incorrect!');
    end
    
end



function c = primal_frobenius(R, t, R1, R2, t1, t2, t_inf, r_inf)

n = size(R1, 3);
d = size(R1, 1);
if nargin < 8 
    t_inf = ones(n,1);
end
if nargin < 7
    r_inf = ones(n,1);
end
c = 0;
for idx=1:n
    c_r = r_inf(idx)*norm(R*R1(:,:,idx) - R2(:,:,idx)*R, 'fro')^2;
    % Make sure the negative signs are right!
    c_t = t_inf(idx)*norm(R*t1(:,idx) + t - (R2(:,:,idx)*t + t2(:,idx)))^2;
    c = c + c_r + c_t;
end

end