function [Ra,Rq, Sk, Ku] = comp_stat(cbv, dn)

Rq = zeros(size(cbv, 1), 1); %store Rq for every cbv
Ra = zeros(size(cbv, 1), 1); %store Ra for every cbv
Sk = zeros(size(cbv, 1), 1); %store Sk for every cbv
Ku = zeros(size(cbv, 1), 1); %store Sk for every cbv

%Evaluate statistics for each codebook vector
for k = 1:size(cbv,1)
    cbv_index = find(dn(:,2)==k);
    dn_cbv = dn(cbv_index, 1); %#ok<FNDSB> 
    Ra(k) = mean(dn_cbv + abs(min(dn_cbv)));
    Rq(k) = rms((dn_cbv + abs(min(dn_cbv))) - Ra(k)); 
%     Rq(k) = sqrt(1/numel(dn_cbv)*sum(dn_cbv.^2));
if abs(Rq(k)) < 2e-6
    Sk(k) = 0;
    Ku(k) = 0;
else
    Sk(k) = skewness(dn_cbv + abs(min(dn_cbv)));
    Ku(k) = kurtosis(dn_cbv + abs(min(dn_cbv)));
    %Sk(k) = (1/(Rq(k)^3))*sum((dn_cbv - Ra(k)).^3);
    %Ku(k) = (1/(Rq(k)^4))*sum((dn_cbv - Ra(k)).^4);
end
end
end