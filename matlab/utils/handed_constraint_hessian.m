function H = handed_constraint_hessian(lam)
%handed_constraint_hessian SO(3) quadratic handedness constraints (i.e. the
% right hand rule) in a hessian for the lagrangian. Derived in Briales
% et al, supplementary material  2017 (CVPR).

n = length(lam);
d = sqrt(n);
assert(d==3, 'Only implemented for SO(3)!');

ids = [1 2 3; 
       2 3 1;
       3 1 2];
   
H = zeros(d^2 + 1, d^2 + 1);

for id=1:d
    i = ids(id, 1);
    j = ids(id, 2);
    k = ids(id, 3);
    eij = zeros(d,d);
    eij(i,j) = 1;
    ek = zeros(d,1);
    ek(k) = 1;
    lam_id = lam((id-1)*d+1:(id-1)*d+3);
    lamskew = skew_matrix(lam_id);
    
    H = H + [-kron(eij, lamskew), -kron(ek, lam_id); 
             zeros(1, 10)];
end

% Make symmetric
H = 0.5*(H + H.');

end

