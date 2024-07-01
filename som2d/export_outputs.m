%% File export routine for the SOM outputs
if wrt_cbv == 1
    writematrix(cbv, append(cbv_filename, '.dat'), 'Delimiter','tab');
end

if wrt_stat == 1
    sb = stat.cbv_arclength;
    ra = stat.mean;
    rq = stat.rms;
    sk = stat.skewness;
    ku = stat.kurtosis;
    tab_stats = table(sb, ra, rq, sk, ku,'VariableNames', {'Arclength', 'Ra', 'Rq', 'Sk', 'Ku'});
    writetable(tab_stats, append(stat_filename, '.dat'), 'Delimiter','tab');
end

if wrt_hmap == 1
    writematrix(hmap.full, append(hmap_filename, '_full.dat'), 'Delimiter','tab');
    if strcmp('stl', hmap_format)
        T = delaunay(hmap.full(:,1),hmap.full(:,3));
        tri = triangulation(T, hmap.full(:,1), hmap.full(:,3), hmap.full(:,2));
        stlwrite(tri, append(hmap_filename, '_full.stl'));
    elseif strcmp('pc', hmap_format)
        pcwrite(pointCloud(hmap.full), append(hmap_filename, '_full.ply'));
    end

elseif wrt_hmap == 2
    writematrix(hmap.ss, append(hmap_filename, '_ss.dat'), 'Delimiter','tab');
    if strcmp('stl', hmap_format)
        T = delaunay(hmap.ss(:,1),hmap.ss(:,3));
        tri = triangulation(T, hmap.ss(:,1), hmap.ss(:,3), hmap.ss(:,2));
        stlwrite(tri, append(hmap_filename, '_ss.stl'));
    elseif strcmp('pc', hmap_format)
        pcwrite(pointCloud(hmap.ss), append(hmap_filename, '_ss.ply'));
    end

elseif wrt_hmap == 3
    writematrix(hmap.ps, append(hmap_filename, '_ps.dat'), 'Delimiter','tab');
    if strcmp('stl', hmap_format)
        T = delaunay(hmap.ps(:,1),hmap.ps(:,3));
        tri = triangulation(T, hmap.ps(:,1), hmap.ps(:,3), hmap.ps(:,2));
        stlwrite(tri, append(hmap_filename, '_ps.stl'));
    elseif strcmp('pc', hmap_format)
        pcwrite(pointCloud(hmap.ps), append(hmap_filename, '_ps.ply'));
    end

elseif wrt_hmap == 4
    writematrix(hmap0.ss, append(hmap0_filename, '_ss.dat'), 'Delimiter','tab');
    writematrix(hmap0.ps, append(hmap0_filename, '_ps.dat'), 'Delimiter','tab');
    if strcmp('stl', hmap0_format)
        T_ss = delaunay(hmap0.ss(:,1),hmap0.ss(:,3));
        tri_ss = triangulation(T_ss, hmap0.ss(:,1), hmap0.ss(:,3), hmap0.ss(:,2));
        stlwrite(tri_ss, append(hmap0_filename, '_ss.stl'));
        T_ps = delaunay(hmap0.ss(:,1),hmap0.ss(:,3));
        tri_ps = triangulation(T_ps, hmap0.ss(:,1), hmap0.ss(:,3), hmap0.ss(:,2));
        stlwrite(tri_ps, append(hmap_filename, '_ps.stl'));
    elseif strcmp('pc', hmap_format)
        pcwrite(pointCloud(hmap0.ss), append(hmap0_filename, '_ss.ply'));
        pcwrite(pointCloud(hmap0.ps), append(hmap0_filename, '_ps.ply'));
    end

end

fprintf('-----> Exported data \n')