function [x, y, z, dz_dx, dz_dy, XX, YY, Z] = random_smooth_traj(N,rad, N_interp, N_modes, ...
                                        amp_range, freq_range)
%random_smooth_traj 

theta = linspace(0, 2*pi, N);
x = cos(theta)*rad; 
y = sin(theta)*rad;

X = linspace(min(x)-1.0, max(x)+1.0, N_interp);
Y = linspace(min(y)-1.0, max(y)+1.0, N_interp);
[XX,YY] = meshgrid(X,Y);
Z = zeros(size(XX));

dz_dx = zeros(size(x));
dz_dy = zeros(size(y));

for idx=1:N_modes
    phase = 2*rand*pi -pi;
    amp = rand*(amp_range(2) - amp_range(1)) + amp_range(1);
    freq = rand*(freq_range(2) - freq_range(1)) + freq_range(1);
    mix = rand;
    Z = Z + amp*cos((mix*XX + (1-mix)*YY)*freq + phase);
    
    dz_dx = dz_dx - amp*sin((mix*x + (1-mix)*y)*freq + phase)*freq*mix;
    dz_dy = dz_dy - amp*sin((mix*x + (1-mix)*y)*freq + phase)*(1-mix)*freq;
    
    % Do it again for Y
%     phase = 2*rand*pi -pi;
%     amp = rand*(amp_range(2) - amp_range(1)) + amp_range(1);
%     freq = rand*(freq_range(2) - freq_range(1)) + freq_range(1);
%     Z = Z + amp*cos(YY*freq + phase);

end



z = interp2(XX,YY,Z,x,y);

% figure;
% surf(XX,YY,Z);

end