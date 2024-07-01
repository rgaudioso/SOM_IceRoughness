function cbv = som_main_val(data, somParams)

tic
%% Create and train the Self-Organizing Map
% SOM Parameters: somParams = {resolution, iter, eta, delta, init}
somSize = somParams.resolution;      % Adjust as needed
epochs = somParams.iter;             % Number of training epochs
eta_initial = somParams.eta;         % Initial learning rate
delta_initial = somParams.delta;     % Initial neighborhood scale parameter
n_weight = somParams.n_weight;             % Weight for linear neighbourhood decreasing function

fprintf('-----> START training loop: Resolution = %d, Epochs = %d \n', somSize, epochs)
cbv = som_train_val(somSize, data, epochs, eta_initial, delta_initial, n_weight);
fprintf('-----> END training loop \n')

% %% Compute metrics
% fprintf('-----> Computing metrics... \n')
% [dn, sb, sx] = comp_norm_arclength(cbv, data);
% [Ra, Rq, Sk, Ku] = comp_stat(cbv, dn);
% sb = sb - sb(find(abs(cbv(:,2))==min(abs(cbv(:,2))))); % translate sb
% stat = struct('cbv_arclength', sb, 'mean', Ra, 'rms', Rq, 'skewness', Sk, 'kurtosis', Ku);

toc

end