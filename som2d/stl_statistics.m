function [hint, stats] = stl_statistics(input_file)
%-------------------------------------------------------------------------%
% This function is intended for post-processing. The input file should be
% already the surface under analysis in terms of dimensions and roughness
% scales to be targeted.
%-------------------------------------------------------------------------%

%% Read inputs !!check units: mm
file = stlread(input_file);
%height map: k = k(x,y)
x = file.Points(:,1);
y = file.Points(:,3);
k = file.Points(:,2);
ka = k + abs(min(k));
% xyk = [x y k];

%% Interpolation of the height map matrix --> try to obtain gridded data for statistics
tol = 0.0001; %set this parameter to a low value (i.e. 1e-4) if the griddata function gives problem with extremal points
dx = max(x)-min(x);
dy = max(y)-min(y);
ds = mean(abs(k))/4; % streamwise and spanwise resolution for interpolation based on mean k (can be adjusted)
nx = floor(dx/ds);
ny = floor(dy/ds);
xq = linspace(min(x)+tol, max(x)-tol, nx); 
yq = linspace(min(y)+tol, max(y)-tol, ny);
[xx, yy] = meshgrid(xq,yq);

% Interpolate height values onto the grid
kk = griddata(x, y, k, xx, yy, 'cubic');
kka = griddata(x, y, ka, xx, yy, 'cubic');
% Plot the interpolated height map
% surf(xx, yy, kk, 'EdgeColor', 'none'), axis equal;
% title('Interpolated Height Map');
% xlabel('X');
% ylabel('Y');
% zlabel('Height (k)');

%% Compute relevant statistics and parameters
% !!check units: mm
% kk = kk.*1000; kka = kka.*1000;
hint = kk;
% Area & hmap dimensions (M,N)
Ap = dx*dy;
M = size(kk,1); N = size(kk,2);
Dx = dx/N; Dz = dy/M;

% Statistics (moments)
ka = (1/(M*N))*sum(abs(kka),'all');
krms = sqrt((1/(M*N))*sum(kk.^2,'all'));
kz = max(kk,[],'all')-min(kk,[],'all');
sk = (1/krms^3)*(1/(M*N))*sum(kk.^3,'all');
ku = (1/krms^4)*(1/(M*N))*sum(kk.^4,'all');

%!!!ACHTUNG: x-wise means in this case along columns, meaning that the corresponding matrix index is the 2nd (j)
%            y-wise means in this case along rows, meaning that the corresponding matrix index is the 1st (i)

% Autocorrelation function using xcorr2
Rk = xcorr2(kk); %(*1/krms^2)-->areal?? See Thakkar et al. 2017 Appendix
Rk_n = Rk/Rk(ceil(end/2), ceil(end/2));

% Visualize the elevation matrix and its autocorrelation
% figure()
% subplot(2, 2, 1);
% imagesc(kka);
% axis equal tight;
% title('Elevation Matrix');
% subplot(2, 2, 2);
% imagesc(Rk_n);
% axis equal tight;
% title('Areal Autocorrelation Matrix');
% subplot(2, 2, 3);
% plot(Rk_n(ceil(end/2), :));
% title('Areal Autocorrelation along X-axis');
% subplot(2, 2, 4);
% plot(Rk_n(:, ceil(end/2)));
% title('Areal Autocorrelation along Y-axis');

% Correlation lengths based on a 0.2 threshold
% lx_corr == minimum lengt at wich the autocorrelation along x drops under 0.2
% lx_corr == minimum lengt at wich the autocorrelation along y drops under 0.2
corr_treshold = 0.2;
i1 = min(find(Rk_n(ceil(end/2),:)>corr_treshold));
i2 = max(find(Rk_n(ceil(end/2),:)>corr_treshold));
j1 = min(find(Rk_n(:,ceil(end/2))>corr_treshold));
j2 = max(find(Rk_n(:,ceil(end/2))>corr_treshold));

if abs(ceil(size(Rk_n,2)/2)-i1) < abs(ceil(size(Rk_n,2)/2)-i2)
    il = i1;
else
    il = i2;
end
if abs(ceil(size(Rk_n,1)/2)-j1) < abs(ceil(size(Rk_n,1)/2)-j2)
    jl = j1;
else
    jl = j2;
end

lx_corr = il*ds; 
ly_corr = jl*ds;
% Surface Anisotropy Ratio --> see Busse et al. 
SAR     = lx_corr/ly_corr;

% Effective slope --> Napoli et al.
[dkx, dky] = gradient(kk,ds,ds);
ESx = (1/(Ap))*trapz(yq,trapz(xq,abs(dkx),2));
ESy = (1/(Ap))*trapz(yq,trapz(xq,abs(dky),2));

% Flack&Schultz ks
kseq = 4.43*krms*(1+sk)^(1.37);

% Max/Min wavelengths
[~, peak_locs_x] = findpeaks(kk(floor(size(kk,1)/2),:));
[~, valley_locs_x] = findpeaks(-kk(floor(size(kk,1)/2),:));
[~, peak_locs_z] = findpeaks(kk(:,floor(size(kk,2)/2)));
[~, valley_locs_z] = findpeaks(-kk(:,floor(size(kk,2)/2)));

wavelength_x = mean(diff(sort([peak_locs_x, valley_locs_x]))).*Dx;
wavelength_z = mean(diff(sort([peak_locs_z', valley_locs_z']))).*Dz;

min_wavelength = min(wavelength_x, wavelength_z);
max_wavelength = max(wavelength_x, wavelength_z);

%% Store results
stats = struct('s', xq, 'ka', ka, 'kz', kz, 'krms', krms, 'Sk', sk, 'Ku', ku, ...
               'Lx_corr', lx_corr, 'Ly_corr', ly_corr, 'SAR', SAR, 'dky',dky,'dkx',dkx,'ESx', ESx, 'ESy', ESy, 'kseq', kseq, ...
                'lambda0', min_wavelength, 'lambda1', max_wavelength);

%% Statistics based on stl area comput.
% [area, area_p] = stl_area(file);
% 
% km = 1/area_p*trapz(yl,trapz(xl,abs(kk),2));
% kt = max(kk,[],'all')-min(kk,[],'all');
% krms = sqrt(1/area_p*trapz(yl,trapz(xl,(kk-km).^2,2)));
% [dkx, dky] = gradient(kk,1/length(xl),1/length(yl));
% ES = 1/area_p*trapz(yl,trapz(xl,abs(dky),2)); %dky or dkx???
% sk = 1/krms^3*trapz(yl,trapz(xl,(kk-km).^3,2));
% ku = 1/krms^4*trapz(yl,trapz(xl,(kk-km).^4,2));

% function to compute the area of the facets
function [area, area_p] = stl_area(file)
area = 0; %initialize the area
area_p = 0;
%loop over all the desired small triangles and accumulate the areas
for i = 1:size(file.ConnectivityList,1)   %the number of rows gives the number of triangles produced
    a =  file.Points(file.ConnectivityList(i,:),:); %this gives the 3 vertices of the ith triangle
    %extract the three vertices
    p1 = a(1,:);  
    p2 = a(2,:);
    p3 = a(3,:);
    %extract the normal
    n = file.faceNormal(i);
    ni = n(1); nj = n(2); 
    area = area + 0.5 * norm(cross(p2-p1,p3-p1));  %calculate the area of the small triangle and add with previous result
    area_p = area_p + 0.5 * norm(cross(p2-p1,p3-p1)) * sin(atan2(ni,nj)); %calculate projected area
end
end

end