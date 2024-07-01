function [init_box, conn] = boxmesh(data, numX, numY, numZ)

nx = numX-1; ny = numY-1; nz = numZ-1;

% Define the range and number of points for each dimension
xmin = min(data(:,1)); xmax = max(data(:,1));
ymin = min(data(:,2)); ymax = max(data(:,2)); 
zmin = min(data(:,3)); zmax = max(data(:,3));
dx = (xmax-xmin)/nx;
dy = (ymax-ymin)/ny; 
dz = (zmax-zmin)/nz;

% Create vectors of evenly??? spaced points
x = linspace(xmin, xmax, numX)';
y = linspace(ymin, ymax, numY)';
z = linspace(zmin, zmax, numZ)';

% box for the airfoil
face1 = []; %xz plane, suction side
for i = 1:nx+1
    xi = xmax - (i-1)*dx; %mind the ordering!!!
    for k = 1:nz+1
        zi = zmin + (k-1)*dz;
        face1 = [face1; xi, ymax, zi];
    end
end
% face1 = [(1:size(face1,1))', face1];

face2 = []; %yz plane, front
for j = 1:ny+1
yi = ymax - (j-1)*dy; %mind the ordering!!!
for k = 1:nz+1
zi = zmin + (k-1)*dz;
face2 = [face2; xmin, yi, zi];
end
end
face2 = face2(numZ+1:end-numZ,:); % exclude repeated edge pts
% face2 = [(size(face1,1)+1:size(face1,1)+size(face2,1))', face2];

% box for the airfoil
face3 = []; %xz plane, pressure side
for i = 1:nx+1
xi = xmin + (i-1)*dx; %mind the ordering!!!
for k = 1:nz+1
zi = zmin + (k-1)*dz;
face3 = [face3; xi, ymin, zi];
end
end
% face3 = [(size(face1,1)+size(face2,1)+1:size(face1,1)+size(face2,1)+size(face3,1))', face3];

init_box = [face1; face2; face3];

num_i = numX + numY + numX -2;
num_j = numZ;
conn  = quad_conn(num_i, num_j);

function conn = quad_conn(num_i, num_j)
ni = num_i-1;
nj = num_j-1;
conn = [];
for ii = 1:ni
    for jj = 1:nj
        id = jj + (ii-1)*nj;
        v1 = id + ii-1;
        v2 = v1 + 1;
        v3 = v2 + num_j;
        v4 = v3 - 1;
        conn = [conn; id v1 v2 v3 v4];
    end
end
end

end