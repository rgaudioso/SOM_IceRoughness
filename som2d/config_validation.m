% !!! Always add the path pointing to the 2d som code !!!
addpath('I:\Documents\TCF\SOM_Roughness\Code\som2d');

%% Sinusoidal roughness test
np = 1e4; % size of the starting analytical dataset
xs = linspace(0,10,np);
ys = 0.25.*sin(0.5*pi.*xs); % sinusoidal profile
%ys = 0.5.*sin(pi.*xs); % sinusoidal profile
L = max(xs)-min(xs);

%% SOM test applied to sinusoidal roughness
data_sin = [xs' ys'];
resolution = 100; %100;
iter = 20;
eta = 0.1;
delta = 1;
n_weight = 0.1;
somParams = struct('resolution', resolution, 'iter', iter, 'eta', eta, 'delta', delta, 'n_weight', n_weight);

%------------------------------RUN -->VALIDATION<-- CASE-----------------------------------%
cbv = som_main_val(data_sin, somParams);

% compute statistics and error wrt analytical profile
Sa = 1/L.*trapz(xs,abs(ys-mean(ys))); % Arithmetic mean height
Ra = 1/resolution*sum(abs(cbv(:,2) - mean(cbv(:,2))));
err = abs((Sa-Ra)/Sa).*100;

% Plot results
figure(), hold on, grid on
plot(xs, ys, 'k-', 'LineWidth', 1.5)
plot(cbv(:,1), cbv(:,2), '-^', 'LineWidth', 1.5)