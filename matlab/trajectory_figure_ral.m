%% Generate data for RAL figure
close all; clear all; clc;
rand('seed', 1);
[r,p,w]=bot_quat_to_roll_pitch_yaw(bot_angle_axis_to_quat(pi/5, [.4 .5 .6]));
trueCalib_euler = [0.6 0.3 0.4 r p w ];
sensor2_expressedIn_sensor1 = GetHomoTransform(trueCalib_euler);


N = 50; % Num poses
N_interp = 300;
N_modes = ceil(N/10);
rad = 10;
amp_range = [0.5, 1];
freq_range = [0.5, 1.5];

[x, y, z,dz_dx,dz_dy, XX, YY, Z] = random_smooth_traj(N,rad, N_interp, N_modes, amp_range, freq_range);
xaxis = [diff(x); diff(y); diff(z)];
yaxis = zeros(size(xaxis));
zaxis = zeros(size(xaxis));
% normal = zeros(size(xaxis));
sensor1_expressedIn_world = zeros(length(x),6);
sensor2_expressedIn_world = zeros(length(x),6);
T = eye(4);
for i = 1:length(x)-1
    T(1:3,1:3) = [xaxis(:,i) yaxis(:,i) zaxis(:,i)]; % * bot_quat_to_matrix(bot_angle_axis_to_quat(pi/3, [1 1 1]));
    T(1:3,4) = [x(i), y(i), z(i)];
    hold on;
    DrawAxis(T, 0.1, 'r', 'b', 'k');
    DrawAxis(T * sensor2_expressedIn_sensor1, 0.1, 'g', 'y', 'c');
    sensor1_expressedIn_world(i,:) = GetState(T);
    sensor2_expressedIn_world(i,:) = GetState(T * sensor2_expressedIn_sensor1);
end
sensor2_expressedIn_world(end,:) = sensor2_expressedIn_world(1,:);

% R_calib = euler2rotm([r,p,w]);
% t_calib = trueCalib_euler(1:3).';
x2 = zeros(size(x));
y2 = zeros(size(y));
z2 = zeros(size(z));
for idx = 1:length(x)
    x2(idx) = sensor2_expressedIn_world(idx,1);
    y2(idx) = sensor2_expressedIn_world(idx,2);
    z2(idx) = sensor2_expressedIn_world(idx,3);
end

xaxis = [diff(x); diff(y); diff(z)];
yaxis = zeros(size(xaxis));
zaxis = zeros(size(xaxis));
normal = zeros(size(xaxis));

%     dz_dx = cos(x);
%     dz_dy = -sin(y);

for i = 1:length(x)-1
    normal(:,i) = cross( [1 0 dz_dx(i)], [0 1 dz_dy(i) ] );
    normal(:,i) = normal(:,i) / norm(normal(:,i));
    xaxis(:,i) = xaxis(:,i) /  norm(xaxis(:,i));        
    yaxis(:,i) = cross ( normal(:,i), xaxis(:,i) );
    yaxis(:,i) = yaxis(:,i) /  norm(yaxis(:,i));
    zaxis(:,i) = cross ( xaxis(:,i), yaxis(:,i) );
    zaxis(:,i) = zaxis(:,i) /  norm(zaxis(:,i));
end

save('../data/trajectory_fig_ral.mat');