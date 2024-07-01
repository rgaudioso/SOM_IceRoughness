function [normalVectors, areas, g] = qnorm_area(vertices, el_index)
    % Calculate normal vectors for quadrilateral surface elements
    numQuads = length(el_index);
    normalVectors = zeros(numQuads, 3);
    areas = zeros(numQuads, 1);
    g = zeros(numQuads,3);

    for i = 1:numQuads
        i0 = 4*(i-1)+1;
        % Choose three non-collinear vertices
        v1 = vertices(i0, :) - vertices(i0+1, :);
        v2 = vertices(i0+3, :) - vertices(i0, :);
        g(i,:) = [mean(vertices(:,1)), mean(vertices(:,2)), mean(vertices(:,3))];
        % Compute the cross product to get the normal vector
        normalVectors(i, :) = cross(v1, v2);
        % Normalize the normal vector
        normalVectors(i, :) = normalVectors(i, :) / norm(normalVectors(i, :));
        % Compute the area as half of the magnitude of the cross product
        areas(i) = 0.5 * norm(normalVectors(:,i));
    end
end