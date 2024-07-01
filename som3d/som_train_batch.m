function cbv = som_train_batch(somSize, input_data, epochs, eta_initial, delta_initial, init_par)

%% CBV initialization !!!
%    % Initialize codebook vectors randomly
%      cbv = rand(somSize, size(input_data, 2));
    
    if nargin==6
        % Initialize codebook vectors with parabola (2d)
        y0 = init_par(1); y1 = init_par(2); a = init_par(3); b = init_par(4); c = init_par(5);
        dy = abs((y1-y0))/(somSize-1);
        yi = (y0:dy:y1)';
        xi = a*(yi-c).^2+b;
        cbv = sortrows([xi, yi],2);
    elseif nargin==5
        cbv = rand(somSize,size(input_data, 2));
    end

    delta = delta_initial;
% Parameters for batch training
numBatches = ceil(size(input_data, 1) / batchSize);
idx = randperm(size(input_data,1)); % Random batches

%% Training loop
for k = 1:numBatches

    batch_idx = idx((k-1)*batchSize + 1 : k*batchSize);
    shuffled_data = input_data(batch_idx, :);

    for epoch = 1:epochs
        
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
        delta = delta * exp(-epoch / epochs);
        % Sort codebook vectors along the Daisy chain
%         cbv = sortrows(cbv,2);
        % Display progress
        fprintf('Trained on epoch %d / %d \n', epoch, epochs)
    end
end

end