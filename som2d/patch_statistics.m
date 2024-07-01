function [hint, stats] = patch_statistics(input_file, npx, npy, res)
%-------------------------------------------------------------------------%
% This function is intended for post-processing. The input file should be
% already the surface under analysis in terms of dimensions and roughness
% scales to be targeted.
%-------------------------------------------------------------------------%
tic
    %% Read inputs !!check units: mm
    file = stlread(input_file);
    % height map: k = k(x,y)
    x = file.Points(:,1);
    y = file.Points(:,3);
    k = file.Points(:,2);
    
    % Determine the range of x and y
    xmin = min(x);
    xmax = max(x);
    ymin = min(y);
    ymax = max(y);

    % Calculate the side length of each subset square
    Lx = (xmax - xmin) / npx;
    Ly = (ymax - ymin) / npy;

    % Initialize statistics structure
    stats = struct();
    figure(), hold on, grid on, axis equal
   
    % Iterate over each subset
    for i = 1:npx
        for j = 1:npy
            % Define subset boundaries
            subset_x_min = xmin + (i - 1) * Lx;
            subset_x_max = min(xmax, xmin + i * Lx);
            subset_y_min = ymin + (j - 1) * Ly;
            subset_y_max = min(ymax, ymin + j * Ly);
            
            % Extract subset of data
            subset_indices = find(x >= subset_x_min & x <= subset_x_max & ...
                                  y >= subset_y_min & y <= subset_y_max);
            subset_x = x(subset_indices);
            subset_y = y(subset_indices);
            subset_k = k(subset_indices);

            %% Interpolation of the height map matrix --> try to obtain gridded data for statistics
            tol = 0.0; % set this parameter for cubic interp to a low value (i.e. 1e-4) if the griddata function gives problem with extremal points
            dx = max(subset_x)-min(subset_x);
            dy = max(subset_y)-min(subset_y);
            ds = mean(abs(subset_k))/res; % streamwise and spanwise resolution for interpolation based on mean k (can be adjusted)
            nx = floor(dx/ds);
            ny = floor(dy/ds);
            xq = linspace(min(subset_x)+tol, max(subset_x)-tol, nx); 
            yq = linspace(min(subset_y)+tol, max(subset_y)-tol, ny);
            [xx, yy] = meshgrid(xq, yq);

            % Interpolate height values onto the grid
            kk = griddata(subset_x, subset_y, subset_k, xx, yy, 'nearest');

            surf(xx,yy,kk,'EdgeColor','none')

            %% Compute relevant statistics and parameters
            % !!check units: mm
            kk = kk.*1e3;
            hint = kk;
            % Area & hmap dimensions (M,N)
            Ap = dx*dy*1e6;
            M = size(kk,1); N = size(kk,2);

            % Statistics (moments)
            ka = (1/(M*N))*sum(abs(kk),'all');
            krms = sqrt((1/(M*N))*sum(kk.^2,'all'));
            kz = max(kk,[],'all')-min(kk,[],'all');
            sk = krms^(-3)*(1/(M*N))*sum(kk.^3,'all');
            ku = krms^(-4)*(1/(M*N))*sum(kk.^4,'all');

            % Autocorrelation function using xcorr2
            Rk = xcorr2(kk);
            Rk_n = Rk/Rk(ceil(end/2), ceil(end/2));

            % Correlation lengths based on a 0.2 threshold
            corr_threshold = 0.2;
            i1 = min(find(Rk_n(ceil(end/2),:) > corr_threshold));
            i2 = max(find(Rk_n(ceil(end/2),:) > corr_threshold));
            j1 = min(find(Rk_n(:,ceil(end/2)) > corr_threshold));
            j2 = max(find(Rk_n(:,ceil(end/2)) > corr_threshold));

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

            lx_corr = il * ds; 
            ly_corr = jl * ds;
            % Surface Anisotropy Ratio --> see Busse et al. 
            SAR = lx_corr / ly_corr;

            % Effective slope --> Napoli et al.
            [dkx, dky] = gradient(kk,ds,ds);
            ESx = (1/(Ap))*trapz(yq,trapz(xq,abs(dkx),2));
            ESy = (1/(Ap))*trapz(yq,trapz(xq,abs(dky),2));

            % Flack&Schultz ks
            kseq = 4.43 * krms * (1 + abs(sk))^(1.37);

            %% Store results for this subset
            subset_stats = struct('ka', ka, 'kz', kz, 'krms', krms, 'Sk', sk, 'Ku', ku, ...
                                  'Lx_corr', lx_corr, 'Ly_corr', ly_corr, 'SAR', SAR, 'ESx', ESx, 'ESy', ESy, 'kseq', kseq);
            
            % Store subset statistics in a cell array
            fieldname = sprintf('subset_%d_%d', i, j);
            stats.(fieldname) = subset_stats;
        end
    end
    toc
end