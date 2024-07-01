function [dn, sb, sx] = comp_norm_arclength(cbv, dataset)

dn = zeros(size(dataset,1), 2);   %store normal distance AND bmu_index
gamma = zeros(size(dataset,1),1); %store local direction of yj points wrt the manifold
sx = zeros(size(dataset, 1), 1);  %store points arclengths
sb = zeros(size(cbv, 1), 1);      %store cbv arclengths

% Evaluate the surface distance coordinate for each codebook vector

sb = zeros(size(cbv, 1), 1);
for k = 2:size(cbv, 1)
    sb(k) = sb(k - 1) + sqrt((cbv(k, 1) - cbv(k - 1, 1))^2 + (cbv(k, 2) - cbv(k - 1, 2))^2);
end

for j = 1:size(dataset, 1)
    xj = dataset(j,:);
    % Find the BMU
    distances = pdist2(xj, cbv);
    [~, bmu_index] = min(distances);
    b = cbv(bmu_index, :);
    
    % Find the neighboring codebook vectors along the "daisy-chain"
    b_minus_1 = cbv(max(1, bmu_index - 1), :);
    b_plus_1 = cbv(min(size(cbv, 1), bmu_index + 1), :);

    % Calculate the direction of the manifold through bn
    alfa_bn = atan2(b_plus_1(2) - b_minus_1(2), b_plus_1(1) - b_minus_1(1));

    % Calculate the direction of the xj point from its winning codebook vector relative to the line through bn
    gamma(j) = atan2(xj(2) - b(2), xj(1) - b(1)) - alfa_bn;

    % Calculate the normal height of the xj point from the line through its winning codebook vector
    dn(j,1) = sqrt((xj(1) - b(1))^2 + (xj(2) - b(2))^2) * sin(gamma(j));
    dn(j,2) = bmu_index;

    % Evaluate the surface projection along the manifold
    sx(j) = sb(bmu_index) + sqrt((xj(1) - b(1))^2 + (xj(2) - b(2))^2) * cos(gamma(j));
end

end