
ni = somSize-1;
nj = n_ds-1;
conn = [];
for jj = 1:nj
    for ii = 1:ni
        id = ii + (jj-1)*ni;
        v1 = id + jj - 1;
        v2 = v1 + ni + 1;
        v3 = v2 + 1;
        v4 = v1 + 1;
        conn = [conn; id v1 v2 v3 v4];
    end
end
n_elems = size(conn,1);
  n_cbv = size(cbv, 1);
    n_pts = size(data, 1);
    
    % Initialize normalDistances
    normalCbv = zeros(n_cbv,3);
    normalElements = zeros(n_elems, 3); baricenter = zeros(n_elems,3); areas = zeros(n_elems,1);
    dn = zeros(n_pts, 1); s = zeros(n_pts, 1);
    weightedNormals = zeros(n_cbv,3);

for i = 1:n_elems
    vertices = cbv(conn(i, 2:end),:);
    % Choose three non-collinear vertices
    v1 = vertices(1, :) - vertices(2, :);
    v2 = vertices(4, :) - vertices(1, :);
    g = [mean(vertices(:,1)), mean(vertices(:,2)), mean(vertices(:,3))];

    % Compute the cross product to get the normal vector
    normalVector = cross(v1, v2);

    % Normalize the normal vector
    normalVector = normalVector/ norm(normalVector);

    % Compute the area as half of the magnitude of the cross product
    area = 0.5 * norm(normalVector);

    % Store results
    areas(i) = area;
    normalElements(i,:) = normalVector;
    baricenter(i,:)    = g;
end

% Iterate over each CBV
 for i_cbv = 1:n_cbv
%       % Find surface elements attached to the current CBV
        [el_index, ~] = find(conn(:, 2:end) == i_cbv);

        % Extract the normal vectors and areas of the attached surface elements
        el_normal = normalElements(el_index,:);
        el_area = areas(el_index,:);

        % Compute the area-weighted mean normal vector
        weightedNormals(i_cbv,:) = mean(el_normal .* el_area);
 end

 for j = 1:n_pts
    xj = data(j,:);
    % Find the BMU
    distances = pdist2(xj, cbv);
    [~, bmu_index] = min(distances);
    
    % Calculate the distance from the point to a plane
    xj_dist = abs(dot(xj - cbv(bmu_index,:), weightedNormals(bmu_index,:)));
%     A = weightedNormals(bmu_index,1); B = weightedNormals(bmu_index,2); C = weightedNormals(bmu_index,3);
%     D = -(cbv(bmu_index,1) + cbv(bmu_index,2) + cbv(bmu_index,3));
%     
%     xj_dist = (abs(A*cbv(bmu_index,1) + B*cbv(bmu_index,2) + C*cbv(bmu_index,3)) +D)/sqrt(A^2 + B^2 + C^2);
    % Evaluate the surface projection along the manifold
%     s_j = sb_vec(bmu_index) + norm(cross(xj - cbv(bmu_index,:), weightedNormals(bmu_index,:)));

    % Store such distance
    dn(j) = xj_dist;
%     s(j)  =s_j;
end
