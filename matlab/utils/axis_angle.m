function R = axis_angle(a, t)
%axis_angle Get a rotation matrix from an axis angle representation
d = length(a);
R = cos(t)*eye(d) + sin(t)*skew_matrix(a) + (1-cos(t))*(a*a.');
end

