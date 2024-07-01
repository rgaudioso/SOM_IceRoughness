function [patch, ss_patch, ps_patch, patch0, ss_patch0, ps_patch0] = shape_unwrap(y, z, sx, dn)

% Find le arc length
s_le = mean(sx(find(abs(y)==min(abs(y)))));

% Translate the arc length 
sx = sx - s_le;

%% Build the patch matrix -> pure roughness, no thickness!
% Suction side
ss_patch = [sx(y>=min(abs(y))) dn(y>=min(abs(y))) z(y>=min(abs(y)))];
% Pressure side
ps_patch = [sx(y<min(abs(y))) dn(y<min(abs(y))) z(y<min(abs(y)))];
% Full patch
patch = [ps_patch; ss_patch];

%% Elevation map -> only roughness peaks ABOVE the manifold
ps_patch0 = ps_patch; ss_patch0 = ss_patch;
for ip = 1:size(ps_patch,1)
    if ps_patch(ip, 2) < 0
        ps_patch0(ip,2) = 0;
    end
end
for is = 1:size(ss_patch,1)
    if ss_patch(is, 2) < 0
        ss_patch0(is,2) = 0;
    end
end
patch0 = [ps_patch0; ss_patch0];

end