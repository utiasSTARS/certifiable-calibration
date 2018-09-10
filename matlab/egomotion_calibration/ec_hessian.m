function H = ec_hessian(muu1,gam,muu2,lam, lam_only)
%ec_hessian 
d = size(muu1,1);

if nargin < 5
    lam_only = false;
end

H_orthog1 = orthog_constraint_hessian(muu1);
H = [H_orthog1 zeros(d^2, 1); zeros(1, d^2) -trace(muu1)-gam];
    

%% If including X.'*X (duality strengthening)
if nargin > 2 && ~lam_only
    H_orthog2 = orthog_constraint_hessian(muu2, false);
    H = H + [H_orthog2 zeros(d^2, 1); zeros(1, d^2) -trace(muu2)];
end

if nargin > 3 
    H_handed = handed_constraint_hessian(lam);
    H = H + H_handed;
end

