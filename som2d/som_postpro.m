function [manifold_plt, ra_plt, rq_plt, sk_plt, ku_plt, hmap_scatter, hmap0_scatter] = som_postpro(geom, cbv, stat, hmap)

%% Interpreter
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

%% Assign data
sb = stat.cbv_arclength;
% sb = (sb+abs(min(sb)))./abs(max(sb)-min(sb)); %only for normalised s between 0 and 1
ra = stat.mean;
rq = stat.rms;
sk = stat.skewness;
ku = stat.kurtosis;

% for roughness map
hmap_full = hmap.full;
% for elevation map
hmap0_full = hmap_full;
for ii = 1:size(hmap_full,1)
    if hmap_full(ii,2)<0
        hmap0_full(ii,2) = 0;
    end
end
% hmap_ss   = hmap.ss;
% hmap_ps   = hmap_ps;

%% Visualize Results
% Manifold plot
manifold_plt = figure(1); hold on, grid on, axis equal
scatter(geom(:,1), geom(:,2), '.')
plot(cbv(:,1), cbv(:,2), '-o', 'linew', 2)
legend('Original Point Cloud', 'Mean Ice Shape')
xlabel('x [m]', 'Interpreter', 'latex'), ylabel('y [m]', 'Interpreter', 'latex')

% Statistics plot
ra_plt = figure(2); hold on, grid on
plot(sb, ra, '-', 'linew', 1.5)
xlabel('s [m]', 'Interpreter', 'latex'), ylabel('Ra [m]', 'Interpreter', 'latex')

rq_plt = figure(3); hold on, grid on
plot(sb, rq, '-', 'linew', 1.5)
xlabel('s [m]', 'Interpreter', 'latex'), ylabel('Rq [m]', 'Interpreter', 'latex')

sk_plt = figure(4); hold on, grid on
plot(sb, sk, '-', 'linew', 1.5)
xlabel('s [m]', 'Interpreter', 'latex'), ylabel('Sk [-]', 'Interpreter', 'latex')

ku_plt = figure(5); hold on, grid on
plot(sb, ku, '-', 'linew', 1.5)
xlabel('s [m]', 'Interpreter', 'latex'), ylabel('Ku [-]', 'Interpreter', 'latex')


% Unwrapped rough patch (full)
hmap_scatter = figure(); hold on, grid on, axis equal
scatter3(hmap_full(:,1), hmap_full(:,3), hmap_full(:,2), 10, hmap_full(:,2), 'filled')
xlabel('s[m]', 'Interpreter', 'latex'), ylabel('z[m]', 'Interpreter', 'latex'), zlabel('dN [m]', 'Interpreter', 'latex')
title('Unwrapped Ice Shape', 'Interpreter', 'latex')

% % Unwrapped elevation patch (full)
% hmap0_scatter = figure(7); hold on, grid on, axis equal
% scatter3(hmap0_full(:,1), hmap0_full(:,3), hmap0_full(:,2), 10, hmap0_full(:,2), 'filled')
% xlabel('s[m]', 'Interpreter', 'latex'), ylabel('z[m]', 'Interpreter', 'latex'), zlabel('dN [m]', 'Interpreter', 'latex')
% title('Unwrapped Ice Shape', 'Interpreter', 'latex')
end

