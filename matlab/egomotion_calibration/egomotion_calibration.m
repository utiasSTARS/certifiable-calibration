function [R,t] = egomotion_calibration(R1, t1, R2, t2)
%ec_calibration Egomotion calibration.

M = ec_get_data_matrix(R1,t1,R2,t2);
% Get Schur complement
Mrr = M(4:end, 4:end);
Mtt = M(1:3, 1:3);
Mrt = M(1:3, 4:end);
Q = Mrr - Mrt.'*inv(Mtt)*Mrt;
% Use maximal linearly indpendent SO(3) constraints
orthog = true;
handed = true;
[~, dual_sol] = ec_dual_solver(Q, orthog, handed);
R = ec_extract_primal(dual_sol);
t = ec_get_t_opt(R, Mtt, Mrt);
end

