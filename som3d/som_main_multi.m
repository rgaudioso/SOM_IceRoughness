clear
close all

%% Load .stl data (X, Y, Z coordinates)
points = stlread("../Case-241-Scan.stl");
x = points.Points(:,1); y = -points.Points(:,2); z = points.Points(:,3); %y needs a '-' sign

%% Limit the data set
dz = max(z)-min(z); zc = (max(z)+min(z))/2; % find the center spanwise coord
dz_lim = 6;                                 % define the total span to be considered
dz_inc = 1;                                 % define the increment-per-dataset (!Must be a sub-multiple of dz_lim!)
zM = zc+dz_lim/2; zm= zc-dz_lim/2;          % define spanwise range
zh = zm+dz_inc/2:dz_inc:zM-dz_inc/2;        % define center spanwise coord for every section

x = x(z(:)>=zm & z(:)<zM); 
y = y(z(:)>=zm & z(:)<zM);
z = z(z(:)>=zm & z(:)<zM);
data = [x y z];

n_ds = dz_lim/dz_inc;
dsets = zeros(size(data,1), 3, n_ds);
for j = 1:dz_lim/dz_inc
    x_ds = x(z(:)>=zm+dz_inc*(j-1) & z(:)<zm+dz_inc*j);
    y_ds = y(z(:)>=zm+dz_inc*(j-1) & z(:)<zm+dz_inc*j);
    z_ds = z(z(:)>=zm+dz_inc*(j-1) & z(:)<zm+dz_inc*j);
    dsets(1:size(x_ds,1),:,j) = [x_ds y_ds z_ds];
end
clear x_ds y_ds z_ds

%% Loop over every dataset
% Parameters
somSize = 80;
epochs = 100;       % Number of training epochs
eta_initial = 0.1;  % Initial learning rate
delta_initial = 3; % Initial neighborhood scale parameter

% Init figures
figure(77), hold on, grid on, axis equal, title('Datasets from original point cloud and cbv'), xlabel('x [m]'), ylabel('y [m]'), zlabel('z [m]')
figure(88), hold on, grid on, axis equal, title('Unwrapped roughness'), xlabel('s [m]'), ylabel('dn [m]'), zlabel('z [m]')

% Store cbv, normals, arc lengths
cbv = zeros(somSize, 3, n_ds);      %init cbv
sb  = zeros(somSize, 1, n_ds);      %init cbv arc lengths
Rq  = zeros(somSize, 1, n_ds);      %init Rq values
dn_mat  = zeros(size(data,1), 2, n_ds); %init normal dist and bmu index
sx_mat  = zeros(size(data,1), 1, n_ds); %init points arc lengths
patch = []; %init unwrapped patch

for i = 1:n_ds
    % Extract and visualize the data set
    dset_i = dsets(:,:,i); dset_i = dset_i(1:max(find(dset_i(:,1)~=0)),:); %throw away redundant zeros!
    figure(77), scatter3(dset_i(:,1), dset_i(:,2), dset_i(:,3), 5, 'filled')

    % Create and train the Self-Organizing Map
%     cbv_i = som(somSize, dset_i(:,1:2));   
    cbv_i = som_train(somSize, dset_i(:,1:2), epochs, eta_initial, delta_initial);

    % Compute metrics
    [dn_i, sb_i, sx_i] = comp_norm_arclength(cbv_i, dset_i(:,1:2));
    Rq_i = comp_rq(cbv_i, dn_i);

%     % unwrapping
%     [patch_i, ~, ~] = shape_unwrap(dset_i(:,2), dset_i(:,3), sx_i, dn_i);
%     figure(88), scatter3(patch_i(:,1), patch_i(:,3), patch_i(:,2), 5, patch_i(:,2), 'filled'), axis equal

    % Assign and store variables
    cbv(:,:,i) = [cbv_i zh(i).*ones(size(cbv_i,1),1)];
    sb(:,:,i)  = sb_i;
    Rq(:,:,i)  = Rq_i;
    dn_mat(1:size(dn_i,1),:,i) = dn_i;
    sx_mat(1:size(sx_i,1),:,i) = sx_i;
%     patch = [patch; patch_i];

    % Display progress
    fprintf('Trained on dataset %d / %d \n', i, n_ds)
end

%% Reshape outputs
cbv = reshape(permute(cbv, [1 3 2]), somSize*n_ds, 3);
sb_vec = reshape(sb, somSize*n_ds, 1);
dn  = []; 
sx  = [];
for k = 1:n_ds
    dn_k = dn_mat(:,:,i); dn_k = dn_k(1:max(find(dn_k(:,1)~=0)),:); %throw away redundant zeros!
    sx_k = sx_mat(:,:,i); sx_k = sx_k(1:max(find(sx_k(:,1)~=0)),:); %throw away redundant zeros!
    dn = [dn; dn_k];
    sx = [sx; sx_k];
end

%% Visualize results
figure(77);
scatter3(cbv(:,1), cbv(:,2), cbv(:,3), 50, 'r', 'filled');

% figure(), scatter3(patch(:,1), patch(:,3), patch(:,2), 5, patch(:,2), 'filled'), axis equal

%% File export
% wrt_input = input('Write output file? [y/n]:','s');
% 
% if wrt_input == 'y'
%     %Write to pointcloud
%     pcwrite(patch,"unwrapped.ply");
% 
%     %Write to stl
%     T = delaunay(patch(:,1),patch(:,3));
%     tri = triangulation(T, patch(:,1), patch(:,3), patch(:,2));
%     stlwrite(tri,"unwrapped.stl");
% else
%     disp('Avoided file export')
% end
