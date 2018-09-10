%% test_random_smooth_traj


N = 100;
N_interp = 300;
N_modes = 5;
rad = 10;
amp_range = [0.5, 1];
freq_range = [0.5, 1.5];

[x, y, z, dz_dx,dz_dy] = random_smooth_traj(N,rad, N_interp, N_modes, amp_range, freq_range);
figure;
plot3(x,y,z);
hold on;
%     quiver3(x,y,z, zaxis(1,:), zaxis(2,:), zaxis(3,:), 0);
%     quiver3(x,y,z, xaxis(1,:), xaxis(2,:), xaxis(3,:), 0);
hold off;
axis equal;
grid on;