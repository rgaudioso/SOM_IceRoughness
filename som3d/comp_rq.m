function Rq = comp_rq(cbv, dn)

Rq = zeros(size(cbv, 1), 1); %store Rq for every cbv

%Evaluate Rq for each codebook vector
for k = 1:size(cbv,1)
    cbv_index = find(dn(:,2)==k);
    dn_cbv = dn(cbv_index,1); %#ok<FNDSB> 
    Rq(k) = rms(dn_cbv);
end

end