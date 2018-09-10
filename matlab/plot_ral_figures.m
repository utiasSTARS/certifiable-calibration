%% Plot figures for RA-L/ICRA submission (Sept. 10, 2018)

%% Trajectory (Figure 3)
close all; clear all; clc;
load('../data/trajectory_fig_ral.mat');

font_size = 14;
blue_color = [0 116 186]/255;
orange_color = [223 80 35]/255;
fig = figure();
set(fig,'defaulttextinterpreter','latex');
plot3(x,y,z+0.1, 'LineWidth', 2.5, 'Color', blue_color);
hold on;
plot3(x2,y2,z2+0.1, 'LineWidth', 2.5, 'Color', orange_color, 'LineStyle', '-.');
s = surf(XX, YY, Z);
lgnd = legend({'Sensor $a$','Sensor $b$', 'Terrain'}, 'Location', 'NorthWest');
set(lgnd, 'Interpreter', 'Latex','FontSize', font_size);
colormap summer
s.EdgeColor = 'none';
s.FaceAlpha = 0.7;
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','FontSize', font_size, 'Interpreter', 'latex');
ylabel('$y$ (m)','FontSize', font_size, 'Interpreter', 'latex');
zlabel('$z$ (m)','FontSize', font_size, 'Interpreter', 'latex');
hold off;
axis equal;
grid on;

%% Create duality gap (Figure 4)
close all; clear all; clc;
load('../data/constraint_comparison_100.mat');

bar_data = [orthog_handed_count(2:end).' handed_count(2:end).' ...
           orthog_count(2:end).' default_count(2:end).'];
max_count = max(orthog_handed_count);
bar_data = 100*bar_data/max_count;

bar_labels = categorical({'$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'});
bar_labels = reordercats(bar_labels,{'$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'});

font_size = 16;
line_width = 2.5;
subplot_margin = 0.12;
subplot_spacing = 0.16;
fig = figure();
set(fig,'defaulttextinterpreter','latex');


subaxis(2,1,1, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
b = bar(bar_labels, bar_data, 'grouped');
hold on;
ax1 = gca;
b(1).FaceColor = [0 116 186]/255;
b(1).FaceAlpha = 0.5;
b(2).FaceColor = [223 80 35]/255;
b(2).FaceAlpha = 0.5;
b(3).FaceColor = [95 129 54]/255;
b(3).FaceAlpha = 0.5;
b(4).FaceAlpha = 0.5;

% % b(3).FaceColor = ;
% % b(4).FaceColor = ;
grid minor
pbaspect([3 1 1])
set(ax1,'FontSize', font_size, 'TickLabelInterpreter','latex');

[lgnd] = legend({'R+C+H', 'R+H', 'R+C', 'R'}, 'Location', 'SouthWest');
set(lgnd, 'Interpreter', 'Latex','FontSize', font_size-4);
set(lgnd.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8([255;255;255;0.9*255]));
ylim([0 110])
% 
% set(ax1,'FontSize', font_size, 'TickLabelInterpreter','latex');
xlabel('Rotation Magnitude (rad)','FontSize', font_size, 'Interpreter', 'latex')
ylabel('\% Optimal','FontSize', font_size, 'Interpreter', 'latex');

% Create duality gap figure with translation
load('../data/constraint_comparison_translation.mat');

bar_data = [squeeze(orthog_handed_count(:,:,2,:)) squeeze(handed_count(:,:,2,:)) ...
            squeeze(orthog_count(:,:,2,:)) squeeze(default_count(:,:,2,:))];
max_count = max(squeeze(orthog_handed_count(:,:,2,:)));
bar_data = 100*bar_data/max_count;

bar_labels = categorical({'$1$', '$10$', '$100$', '$1000$'});
bar_labels = reordercats(bar_labels,{'$1$', '$10$', '$100$', '$1000$'});

font_size = 16;
line_width = 2.5;
% fig = figure();
% set(fig,'defaulttextinterpreter','latex');
subaxis(2,1,2, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
b = bar(bar_labels, bar_data, 'grouped');
hold on;
ax1 = gca;
b(1).FaceColor = [0 116 186]/255;
b(1).FaceAlpha = 0.5;
b(2).FaceColor = [223 80 35]/255;
b(2).FaceAlpha = 0.5;
b(3).FaceColor = [95 129 54]/255;
b(3).FaceAlpha = 0.5;
% b(4).FaceColor = ;
b(4).FaceAlpha = 0.5;

grid minor
pbaspect([3 1 1])
set(ax1,'FontSize', font_size, 'TickLabelInterpreter','latex');

% [lgnd] = legend({'R + C + H', 'R + H', 'R + C', 'R'}, 'Location', 'SouthWest');
% set(lgnd, 'Interpreter', 'Latex','FontSize', font_size-4);
% set(lgnd.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8([255;255;255;0.9*255]));
ylim([0 110])
% 
% set(ax1,'FontSize', font_size, 'TickLabelInterpreter','latex');
xlabel('Translation Magnitude (m)','FontSize', font_size, 'Interpreter', 'latex')
ylabel('\% Optimal','FontSize', font_size, 'Interpreter', 'latex');


%% Create histograms (Figure 5)
close all; clear all; clc;

load('../data/random_smooth_noise_realistic_histogram.mat');
font_size = 14;
num_bins = 15;
subplot_margin = 0.1;
subplot_spacing = 0.02;

blue_color = [0 116 186]/255;
orange_color = [223 80 35]/255;

x0 = 0;
y0 = 0;
width = 8;
height = 10;
fig = figure('Units','inches',...
'Position',[x0 y0 width height]);

for i = 1:2:24
    j = i;
%     if i > 8
%         j = j + 8;
%     end
    if mod(j+1,4) == 0
        rot_index = 2;
    else
        rot_index = 1;
    end

    subaxis(6,4,i, 'Margin', subplot_margin, 'Spacing', subplot_spacing)
    
    
    histogram(hist_tran_err(j+1,:),linspace(0,4,num_bins), 'Normalization', 'Probability', ...
    'EdgeColor', blue_color, 'LineWidth', 0.5, 'FaceColor', blue_color, 'FaceAlpha', 0.2)

    hold on;
    
    histogram(hist_tran_err(j,:),linspace(0,4,num_bins), 'Normalization', 'Probability', ...
        'EdgeColor', orange_color, 'LineWidth', 0.5, 'FaceColor', orange_color, 'FaceAlpha', 0.6)%    [f, x] = histcounts(hist_tran_err(j,:),linspace(0,5,num_bins), 'Normalization', 'cdf');

    hold off;
    grid minor
    xlim([0 4])
    pbaspect([3 2 1])
    set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');
    if i > 20
        xlabel('Trans. Error (m)','FontSize', font_size, 'Interpreter', 'latex')
        xticks(1:4)
    else
        xticks(1:4)
        xticklabels({})
    end
    
    if mod(i, 4) == 1
        %ylabel('Frequency','FontSize', font_size, 'Interpreter', 'latex');
    else
        yticklabels({})
    end
    
    if i == 1
        lgnd = legend({'Convex','Local'}, 'Location', 'NorthEast');
        set(lgnd, 'Interpreter', 'Latex','FontSize', font_size-2);
    end

    
    subaxis(6,4,i+1, 'MT', subplot_margin, 'MB', subplot_margin,'ML', subplot_margin-0.0075, 'MR', subplot_margin+0.0075, 'Spacing', subplot_spacing)

    histogram(hist_rot_err(j+1,:),linspace(0,0.8,num_bins), 'Normalization', 'Probability', ...
    'EdgeColor', blue_color, 'LineWidth', 0.5, 'FaceColor', blue_color, 'FaceAlpha', 0.2);
    hold on;
    
    histogram(hist_rot_err(j,:),linspace(0,0.8,num_bins), 'Normalization', 'Probability', ...
        'EdgeColor', orange_color, 'LineWidth', 0.5, 'FaceColor', orange_color, 'FaceAlpha', 0.6)%    [f, x] = histcounts(hist_tran_err(j,:),linspace(0,5,num_bins), 'Normalization', 'cdf');

    str = ['$$\sigma_\mathrm{t} = ',num2str(iter_tra(ceil(j/4))), '$$', newline,'$$\sigma_\mathrm{r} =',num2str(iter_rot(rot_index)), '$$'];
%     disp(str)
    text(0.37,0.7,str,'Interpreter','latex', 'FontSize', font_size-2, 'BackgroundColor', 'w', 'EdgeColor', 0.6*ones(1,3));
    
    hold off;
    grid minor
    
    xlim([0 0.8])
    yticklabels({})
    
    ylim([0 1])
    
    
    set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');
    pbaspect([3 2 1])
    
    if i > 20
        xlabel(['Rot. Error'],'FontSize', font_size, 'Interpreter', 'latex');
        xticks([0.2 0.4 0.6 0.8]);
        xticklabels({'$0.2$', '$0.4$', '$0.6$', '$0.8$'});
    else
%         xticks([pi/16 pi/8 3*pi/16 pi/4])
        xticks([0.2 0.4 0.6 0.8]);
        xticklabels({})
    end
end

%% Create heatmaps (Figure 6)
close all; clear all; clc;

font_size = 14;
x0 = 0;
y0 = 0;
width = 8;
height = 4;
blue_color = [0 116 186]/255;
orange_color = [223 80 35]/255;
cmap = [blue_color; orange_color];

fig = figure('Units','inches',...
'Position',[x0 y0 width height]);

set(fig,'defaulttextinterpreter','latex');

load('../data/heatmap_high_res_sept3.mat');
subaxis(1,2,1, 'Margin', 0.1, 'Spacing', 0.02)

data = max(err_b_rot_heat,[],3);
%// Define integer grid of coordinates for the above data
[X,Y] = meshgrid(1:size(data,2), 1:size(data,1));

%// Define a finer grid of points
[X2,Y2] = meshgrid(1:0.1:size(data,2), 1:0.1:size(data,1));

%// Interpolate the data and show the output
outData = interp2(X, Y, data, X2, Y2, 'linear');
imagesc(outData);
pbaspect([1,1,1])

colormap('parula')
c = colorbar;
set(c,'FontSize', font_size, 'TickLabelInterpreter','latex');
set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');

ylabel('Intial Trans. Norm (m)', 'FontSize', font_size-2)
xlabel('Intial Rot. Angle', 'FontSize', font_size-2)

x_tick_ids = [1:size(outData, 1)/4:size(outData, 1) size(outData, 1)];
yticks(x_tick_ids)
yticklabels({'0.2','2.7','5.2','7.7','10.2'})


x_tick_ids = [1:size(outData, 1)/4:size(outData, 1) size(outData, 1)];
xticks(x_tick_ids)
xticklabels({'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$'})
title('Rotation Error','FontSize', font_size, 'Interpreter', 'latex');
subaxis(1,2,2, 'Margin', 0.1, 'Spacing', 0.02)

data = max(err_b_tra_heat,[],3);
%// Define integer grid of coordinates for the above data
[X,Y] = meshgrid(1:size(data,2), 1:size(data,1));

%// Define a finer grid of points
[X2,Y2] = meshgrid(1:0.1:size(data,2), 1:0.1:size(data,1));

%// Interpolate the data and show the output
max_alpha = max(alphas);
max_trans = max(trans_mag);

outData = interp2(X, Y, data, X2, Y2, 'linear');
imagesc(outData);
set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');

yticklabels('')

x_tick_ids = [1:size(outData, 1)/4:size(outData, 1) size(outData, 1)];
xticks(x_tick_ids)

xticklabels({'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$'})

xlabel('Intial Rot. Angle', 'FontSize', font_size-2)
title('Translation Error (m)','FontSize', font_size, 'Interpreter', 'latex');
pbaspect([1,1,1])

%// Add colour bar
c = colorbar;
set(c,'FontSize', font_size, 'TickLabelInterpreter','latex');

%% Create timing figure (Figure 7)
close all; clear all; clc;

load('../data/time_vs_N');

font_size = 14;
line_width = 1.5;
blue_color = [0 116 186]/255;
orange_color = [223 80 35]/255;

fig = figure(); 
set(fig,'defaulttextinterpreter','latex');

errorbar(sim_N_iter,mean(solving_times_vsN_b,2),mean(solving_times_vsN_b,2)...
    -quantile(solving_times_vsN_b,0.25,2),...
    quantile(solving_times_vsN_b,0.75,2)-mean(solving_times_vsN_b,2), '-s', ...
'LineWidth', 0.75,'MarkerSize',5,...
    'MarkerEdgeColor',blue_color,'MarkerFaceColor',blue_color)

hold on

errorbar(sim_N_iter,mean(solving_times_vsN_c,2),mean(solving_times_vsN_c,2)...
    -quantile(solving_times_vsN_c,0.25,2),...
    quantile(solving_times_vsN_c,0.75,2)-mean(solving_times_vsN_c,2), '-s', ...
'LineWidth', 0.75,'MarkerSize',5,...
    'MarkerEdgeColor',orange_color,'MarkerFaceColor',orange_color)

plot(sim_N_iter,mean(solving_times_vsN_b,2), 'Color', blue_color,'LineWidth', line_width)
plot(sim_N_iter,mean(solving_times_vsN_c,2), 'Color', orange_color, 'LineWidth', line_width)

hold off;

grid minor
% pbaspect([2 1 1])
set(gca,'FontSize', font_size, 'TickLabelInterpreter','latex');
xlabel('Number of Measurements','FontSize', font_size, 'Interpreter', 'latex')
ylabel('Time (s)','FontSize', font_size, 'Interpreter', 'latex');
pbaspect([3 1 1])
lgnd = legend({'Local','Convex'}, 'Location', 'NorthWest');
set(lgnd, 'Interpreter', 'Latex','FontSize', font_size);
