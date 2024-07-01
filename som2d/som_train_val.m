function cbv = som_train_val(somSize, input_data, epochs, eta_initial, delta_initial, n_weight)


%% CBV initialization !!!
% Random initialization
cbv = [(linspace(min(input_data(:,1)), max(input_data(:,1)), somSize))' zeros(somSize, 1)];
% cbv = sortrows(rand(somSize,2),1);
% Initial neighbourhood size
delta = delta_initial;

%% Training loop
    for epoch = 1:epochs
        % Shuffle data for each epoch
        shuffled_data = input_data(randperm(size(input_data, 1)), :);
        
        % Update learning rate for the current epoch
        eta = eta_initial * exp(-epoch / epochs);
        
        for i = 1:size(shuffled_data, 1)
            x = shuffled_data(i, :);
            
            % Find the winning codebook vector
            [~, win_idx] = min(sum((cbv - x).^2, 2));
            
            % Update codebook vectors based on the neighborhood function
            for j = 1:size(cbv,1)
                h = exp(-(win_idx - j)^2 / (delta^2));
                cbv(j, :) = cbv(j, :) + h * eta * (x - cbv(j, :));
            end
        end
        
        % Update neighborhood scale for the next epoch
%         delta = delta * exp(-epoch / epochs);
        delta = delta * (1 - n_weight*(epoch/epochs));
        % Display progress
        fprintf('Trained on epoch %d / %d \n', epoch, epochs)
    end
end