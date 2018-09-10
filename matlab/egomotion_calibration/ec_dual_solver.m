function [dual_val, dual_sol] = ec_dual_solver(Q,orthog_str, handed_str)
%ec_dual_solver 

if nargin < 2
    orthog_str = false;
end
if nargin < 3
    handed_str = false;
end

d = sqrt(size(Q, 1) - 1);
assert(d == 3, 'Only support SO(3) for now.');

cvx_begin quiet
    variable muu1(d,d) symmetric
    variable muu2(d,d) symmetric
    variable lam(d^2) 
    variable gam

    maximize(gam);
    subject to
         if handed_str && orthog_str
            Q + ec_hessian(muu1,gam,muu2,lam) == semidefinite(d^2 + 1);
         elseif orthog_str
            Q + ec_hessian(muu1,gam,muu2) == semidefinite(d^2 + 1);
         elseif handed_str
            Q + ec_hessian(muu1,gam,muu2,lam, true) == semidefinite(d^2 + 1);
         else
            Q + ec_hessian(muu1,gam) == semidefinite(d^2 + 1);
         end
cvx_end

dual_val = gam;
dual_sol.gam = gam;
dual_sol.muu1 = muu1;

if handed_str
    dual_sol.lam = lam;
end
if orthog_str
    dual_sol.muu2 = muu2;
end

if handed_str
    dual_sol.H = Q + ec_hessian(muu1,gam,muu2,lam);
elseif orthog_str
    dual_sol.H = Q + ec_hessian(muu1,gam,muu2);
else
    dual_sol.H = Q + ec_hessian(muu1,gam);
end

end

