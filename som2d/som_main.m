function [geom, cbv, stat, hmap, hmap0] = som_main(scan, datatype, units, span_lim, somParams)

tic 

%% Load data in .stl or .dat format (datatype 1 or 2)
if datatype==1
    points = stlread(scan);
    x = points.Points(:,1); y = points.Points(:,2); z = points.Points(:,3);
elseif datatype==2
    points = readtable(scan, 'NumHeaderLines', 1);
    x = points.Var1; y = points.Var2; z = points.Var3;
end

% Check the units and convert to m if needed
if ~strcmp('m', units)   
    if strcmp('mm', units)
        mmtom = 0.001;
        x = x.*mmtom; y = y.*mmtom; z = z.*mmtom;
    end
    if strcmp('in', units)
        intom = 0.0254;
%         x = x.*intom; y = -y.*intom; z = z.*intom; % Generally, for IPW cases -y is required
        x = x.*intom; y = y.*intom; z = z.*intom;
    end
end

% Limit the data set if the limit is specified
if span_lim > 0
    %Zc = (max(z) + min(z))/2; dZ_lim = span_lim; %specified span_lim
    Zc = (max(z) + min(z))/2; dZ_lim = (span_lim/100)*(abs(max(z) - min(z))); %perc span_lim
    U = Zc + dZ_lim/2; L = Zc - dZ_lim/2; % Upper and lower spanwise limits
    x = x(z(:)>=L & z(:)<U);
    y = y(z(:)>=L & z(:)<U);
    z = z(z(:)>=L & z(:)<U);
end

data = [x y]; % Matrix of data points, 2d
geom = data;

fprintf('-----> Loaded geometry points \n')

%% Create and train the Self-Organizing Map
% SOM Parameters: somParams = {resolution, iter, eta, delta, init}
somSize = somParams.resolution;      % Adjust as needed
epochs = somParams.iter;             % Number of training epochs
eta_initial = somParams.eta;         % Initial learning rate
delta_initial = somParams.delta;     % Initial neighborhood scale parameter
n_weight = somParams.n_weight;             % Weight for linear neighbourhood decreasing function
init_par = somParams.init;           % (y0,y1,a,b,c, x = a*(y-c)^2+b, y = y0:dy:y1, dy = (y1-y0)/somSize

fprintf('-----> START training loop: Resolution = %d, Epochs = %d \n', somSize, epochs)
cbv = som_train(somSize, data, epochs, eta_initial, delta_initial, n_weight, init_par);
fprintf('-----> END training loop \n')

%% Compute metrics
fprintf('-----> Computing metrics... \n')
[dn, sb, sx] = comp_norm_arclength(cbv, data);
[Ra, Rq, Sk, Ku] = comp_stat(cbv, dn);
if ~isempty(init_par) % I want to know where the le is only for airfoils
sb = sb- sb(find(abs(cbv(:,2))==min(abs(cbv(:,2))))); % translate sb
end
stat = struct('cbv_arclength', sb, 'mean', Ra, 'rms', Rq, 'skewness', Sk, 'kurtosis', Ku);

%% Unwrapping
fprintf('-----> Unwrapping... \n')
[patch, ss_patch, ps_patch, patch0, ss_patch0, ps_patch0] = shape_unwrap(y, z, sx, dn);
hmap = struct('full', patch, 'ss', ss_patch, 'ps', ps_patch);
hmap0 = struct('full', patch0, 'ss', ss_patch0, 'ps', ps_patch0);

toc
end