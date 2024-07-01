function cbv = som_train(somSize, input_data, epochs, eta_initial, delta_initial, n_weight, init_par)

%% CBV initialization !!!

% Initialize codebook vectors randomly
if isempty(init_par)
    cbv = [(linspace(min(input_data(:,1)), max(input_data(:,1)), somSize))' zeros(somSize, 1)];
%     cbv = flip(cbv);
else
    % Initialize codebook vectors with parabola (2d)
    y0 = init_par(1); y1 = init_par(2); a = init_par(3); b = init_par(4); c = init_par(5);
    dy = abs((y1-y0))/(somSize-1);
    yi = (y0:dy:y1)';
    xi = a*(yi-c).^2+b;
    cbv = sortrows([xi, yi],2);
end

% % Initialize codebook vectors with box
% xm = min(input_data(:,1));  xM_up = max(input_data(input_data(:,2)>=0,1)); xM_bot = max(input_data(input_data(:,2)<0,1)); 
% ym = min(input_data(:,2));  yM = max(input_data(:,2));
% dist = abs(xM_up-xm) + abs(xM_bot-xm) + abs(yM-ym);
% ds = dist/somSize;
% cbv_xup = [(xm:ds:xM_up)' yM.*ones(length(xm:ds:xM_up),1)];
% cbv_y = [xm.*ones(length(ym:ds:yM),1) (ym:ds:yM)'];
% cbv_xbot = [flip((xm:ds:xM_bot)') ym.*ones(length(xm:ds:xM_bot),1)];
% cbv = [cbv_xbot; cbv_y(2:end-1,:); cbv_xup];

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