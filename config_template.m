%% Configuration file for SOM clustering

% !!! Always add the path pointing to the 2d som code !!!
addpath('I:\Documents\TCF\SOM_Roughness\Code\som2d');

%------------------------------INPUT SECTION------------------------------%
% Define the inputs for the som_main function. Inputs and their definition
% are briefly described below. The following inputs are required:
% [cbv, stat, hmap] = som_main(scan, datatype, units, span_lim, somParams) 
% - scan --> type: string
%            role: defines the inital geometry file path
% - datatype --> type: float
%                role: specifies the format of the initial geometry file
%                (1 = .stl, 2 = table format (.dat, .csv, .txt ...)
% - units --> type: string
%             role: defines the inital geometry meas. units ('m' = meters,
%             'in' = inches)
% - span_lim --> type: float
%                role: defines the spanwise extension desired for the
%                geometry. If 0, no limitations are applied to the data in
%                spanwise direction. Check the spanwise geometry extension
%                before specifying this parameter
% - somParams --> type: struct
%                 role: structure allocating all parameters required by the
%                 training (som_train). Specifically:
%                 somParams = {resolution = n. of SOM points or cbv, float
%                              iter       = n. of training epochs, float
%                              eta        = init learning rate, float
%                              delta      = init neighbour param, float
%                              n_weight   = weigth for neighbourhood function decrease, float
%                              init       = cbv initialization (quadratic)
%                                           double with entries [y0, y1, a, b, c]
%                                           x = a*(y - c).^2 + b}
%                 !!! init is case-sensitive: adjust the parameters of the
%                 parabola to fit your input data !!!

scan = 'path-to-geometryfile';
datatype = 1;
units = 'mm';
span_lim = 0;
resolution = 80;
iter = 20;
eta = 0.1;
delta = 1;
n_weight = 0.2;
init = [-0.030, 0.050, 50.0, -0.0035, 0.005]; 
somParams = struct('resolution', resolution, 'iter', iter, 'eta', eta, 'delta', delta, 'n_weight', n_weight, 'init', init);

%------------------------------RUN CASE-----------------------------------%
% Run the SOM clustering based on input data specified above and store the
% outputs: 
% [cbv, stat, hmap] = som_main(scan, datatype, units, span_lim, somParams)
% - geom --> type: double
%            role: matrix storing original geometry points [xp, yp, zp]
% - cbv --> type: double
%           role: codebook vectors matrix [x, y]
% - stat --> type: struct
%            role: statistics struct {sb, Ra, Rq, Sk, Ku}
% - hmap --> type: struct
%            role: height map struct {full, ss, ps}; 
%            each rough patch is stored as [s, dn, z]

[geom, cbv, stat, hmap] = som_main(scan, datatype, units, span_lim, somParams);

%-------------------------------EXPORT------------------------------------%
% Specify the outputs to be written in .dat format. The flag 1 writes the 
% output, the flag 0 avoids the writing. The hmap flag 1 writes the full 
% patch, 2 writes only the suction side patch, the flag 3 only the pressure
% side, the flag 4 stores both separately. The output can be written to an
% .stl file representing the unwrapped surface, or to a .ply file storing
% the unwrapped surface cloud. WARNING: don't trust Matlab triangulation...

% Flags and filenames w/o extension!
wrt_cbv  = 1; 
cbv_filename  = 'manifold';
wrt_stat = 1; 
stat_filename = 'statistics';
wrt_hmap = 1; 
hmap_filename = 'hmap'; 
hmap_format = 'stl'; 

% Run the routine to export the specified files
export_outputs;

%--------------------------- POST-PROCESSING------------------------------%
% Post-process the outputs to produce plots, only if the pospro flag is
% set to 1. The plots must be edited and stored outside of this config
% file, which uses standard options. Access the som_pospro function to 
% suppress undesired outputs or modify the appearence for each plot:
% [manifold_plt, ra_plt, rq_plt, sk_plt, ku_plt, hmap_scatter] = som_postpro(geom, cbv, stat, hmap)
% - manifold_plt --> type: figure
%                    role: plot of original pc and cbv manifold
% - ra_plt --> type: figure
%              role: plot of mean roughness vs arc length
% - rq_plt --> type: figure
%              role: plot of rms roughness vs arc length
% - sk_plt --> type: figure
%              role: plot of roughness sknewness vs arc length
% - ku_plt --> type: figure
%              role: plot of roughness kurtosis vs arc length
% - hmap_scatter --> type: figure
%                    role: scatter plot of the unwrapped ice roughness

postpro = 1;

if postpro == 1
    [manifold_plt, ra_plt, ~, ~, ~, hmap_scatter] = som_postpro(geom, cbv, stat, hmap);
end

%-----------------------END of configuration file-------------------------%