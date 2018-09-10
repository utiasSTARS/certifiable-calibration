function [A, b, c] = ec_get_affine_cost_matrices(R1, t1, R2, t2, r_inf, t_inf)
%ec_get_affine_cost_matrices Ri, ti are stacked data matrices of egomotion
% (incremental motions in SE(3)). Returns affine QCQP cost function 
% matrices: f(x) = x.'*A*x + b.'*x + c.
% R1, R2  - (d x d x n) tensor
% t1, t2  - (d x n) matrix

n = size(R1, 3);
d = size(R1, 1);

if nargin < 6 
    t_inf = ones(n,1);
end
if nargin < 5
    r_inf = ones(n,1);
end

A = zeros(d^2 + d, d^2 + d);
b = zeros(d^2 + d, 1);
c = 0;
for idx=1:n
    R1i = R1(:,:,idx);
    R2i = R2(:,:,idx);
    t1i = t1(:, idx);
    t2i = t2(:, idx);
    r_inf_i = r_inf(idx);
    t_inf_i = t_inf(idx);
    Mtt = get_mtt(R2i, t_inf_i);
    Mrr = get_mrr(R1i, R2i, t1i, t2i, t_inf_i, r_inf_i);
    Mrt = get_mrt(R2i, t1i, t2i, t_inf_i);
    A = A + [Mtt Mrt; Mrt.' Mrr];
    b = b + 2*get_b(R1i, R2i, t1i, t2i, t_inf_i, r_inf_i);
    c = c + t2i.'*t2i;
end

end



function Mtt = get_mtt(R2, t_inf)
d = size(R2, 1);
Mtt = t_inf*(eye(d) - R2).'*(eye(d) - R2);
end

function Mrr = get_mrr(R1, R2, t1, t2, t_inf, r_inf)
d = size(R1, 1);
A = kron(R1.', eye(d)) - kron(eye(d), R2);
Mrr = r_inf*(A.'*A) + t_inf*kron(t1.', eye(d)).'*(kron(t1.', eye(d)));
% Mrc = -t_inf*kron(t1, eye(d))*t2;
% Mrr = [Mrr Mrc; Mrc.' t2.'*t2];

end

function Mrt = get_mrt(R2, t1, t2, t_inf)
d = size(R2, 1);
Mrt = (eye(d) - R2).'*kron(t1.', eye(d));
% Mrt = t_inf*[Mrt  (eye(d)- R2).'*(-t2)];
Mrt = t_inf*Mrt;
end

function b = get_b(R1, R2, t1, t2, t_inf, r_inf)
d = size(R1, 1);
Mrc = -t_inf*kron(t1, eye(d))*t2;
Mtc = t_inf*(eye(d)- R2).'*(-t2);

b = [Mtc; Mrc];

end