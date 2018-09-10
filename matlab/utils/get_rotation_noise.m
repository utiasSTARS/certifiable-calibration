function R = get_rotation_noise(sigma)
%get_rotation_noise Returns rotation matrix with random axis and angle with
% normal distribution with st. dev. sigma and zero mean.
a = rand(3,1)*2 - 1;
a = a./norm(a);
R = axis_angle(a, randn(1)*sigma);
end

