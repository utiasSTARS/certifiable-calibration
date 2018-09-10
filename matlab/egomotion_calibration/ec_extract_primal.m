function R = ec_extract_primal(dual_sol)
%ec_extract_primal ADAPTED FROM JESUS BRIALES' code: 
%                  https://github.com/jbriales/CVPR17

[U,S,~] = svd(dual_sol.H);
s = diag(S);

% Check that smallest eig is small enough (i.e. numerical error, should be
% zero)
zero_check = s(end)<1e-3 && (s(end)/s(end-1))<1e-3;

if zero_check
  % Recovering the solution is trivial, just scale the eigenvector
  r = makenonhom(U(:,end));
  R = reshape(r,3,3);
else
  r = makenonhom(U(:,end));
  R = reshape(r,3,3);
  warning('No solution for non-tight case through nullspace yet')
end

end

function nonhom_x = makenonhom(x)
% nonhom_x = makenonhom( x )
% Returns the homogeneous vector (append 1 at the end of a column vector)

if size(x,2) ~= 1
  error('Use only column vectors with makenonhom');
end

if abs(x(end)) < 1e-6
  warning('Hom component is zero');
end

nonhom_x = x(1:end-1) / x(end);

end