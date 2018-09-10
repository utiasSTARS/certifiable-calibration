function H = orthog_constraint_hessian(muu, cols)
%orthog_constraint_hessian 

if nargin < 2
    cols = true;
end

if cols
    H = 0.5*kron(muu, eye(size(muu, 1)));
else
    H = 0.5*kron(eye(size(muu,1)), muu);
end
% Only the off-diagonal needs to be halved - I don't think this affects any
% of the optimization problems, but need to TEST IT!
H = H + diag(diag(H));
end

