function results = run_duality_experiment(R1,t1,R2,t2)
%run_duality_experiment 

M = ec_get_data_matrix(R1,t1,R2,t2);
% Get Schur complement
Mrr = M(4:end, 4:end);
Mtt = M(1:3, 1:3);
Mrt = M(1:3, 4:end);
Q = Mrr - Mrt.'*inv(Mtt)*Mrt;

% Default solution
% disp('Default');
[dual_val, dual_sol] = ec_dual_solver(Q, false, false);
R_dual = ec_extract_primal(dual_sol);
t_dual = ec_get_t_opt(R_dual, Mtt, Mrt);
primal_val = [vec(R_dual); 1].'*Q*[vec(R_dual); 1];
default.dual_val = dual_val;
default.dual_sol = dual_sol;
default.primal_val = primal_val;
default.gap = primal_val - dual_val;
default.R = R_dual;
default.t = t_dual;
results.default = default;

% Orthog only
% disp('Orthog');
[dual_val, dual_sol] = ec_dual_solver(Q, true, false);
R_dual = ec_extract_primal(dual_sol);
t_dual = ec_get_t_opt(R_dual, Mtt, Mrt);
primal_val = [vec(R_dual); 1].'*Q*[vec(R_dual); 1];
orthog.dual_val = dual_val;
orthog.dual_sol = dual_sol;
orthog.primal_val = primal_val;
orthog.gap = primal_val - dual_val;
orthog.R = R_dual;
orthog.t = t_dual;
results.orthog = orthog;

% Handed only
% disp('Handed');
[dual_val, dual_sol] = ec_dual_solver(Q, false, true);
R_dual = ec_extract_primal(dual_sol);
t_dual = ec_get_t_opt(R_dual, Mtt, Mrt);
primal_val = [vec(R_dual); 1].'*Q*[vec(R_dual); 1];
handed.dual_val = dual_val;
handed.dual_sol = dual_sol;
handed.primal_val = primal_val;
handed.gap = primal_val - dual_val;
handed.R = R_dual;
handed.t = t_dual;
results.handed = handed;

% Both
% disp('Orthog Handed');
[dual_val, dual_sol] = ec_dual_solver(Q, true, true);
R_dual = ec_extract_primal(dual_sol);
t_dual = ec_get_t_opt(R_dual, Mtt, Mrt);
primal_val = [vec(R_dual); 1].'*Q*[vec(R_dual); 1];
orthog_handed.dual_val = dual_val;
orthog_handed.dual_sol = dual_sol;
orthog_handed.primal_val = primal_val;
orthog_handed.gap = primal_val - dual_val;
orthog_handed.R = R_dual;
orthog_handed.t = t_dual;
results.orthog_handed = orthog_handed;


end
