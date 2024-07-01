function szn = arclength(dataset)

dataset_ss = sortrows(dataset(dataset(:,2)>=0,:), 1);
dataset_ps = sortrows(dataset(dataset(:,2)<0,:), 1);

s_ss = zeros(size(dataset_ss, 1), 3);
s_ps = zeros(size(dataset_ps, 1), 3);

s_ss(1,1) = sqrt(dataset_ss(1, 1)^2 + dataset_ss(1, 2)^2); s_ss(1, 2) = dataset_ss(1, 3); s_ss(1, 3) = dataset_ss(1, 4);
s_ps(1,1) = -sqrt(dataset_ps(1, 1)^2 + dataset_ps(1, 2)^2); s_ps(1, 2) = dataset_ps(1, 3); s_ss(1, 3) = dataset_ss(1, 4);

for is = 2:size(dataset_ss, 1)
    s_ss(is,1) = s_ss(is - 1) + sqrt((dataset_ss(is, 1) - dataset_ss(is - 1, 1))^2 + (dataset_ss(is, 2) - dataset_ss(is - 1, 2))^2);
    s_ss(is,2) = dataset_ss(is,3);
    s_ss(is,3) = dataset_ss(is,4);
end
for ip = 2:size(dataset_ps, 1)
    s_ps(ip,1) = s_ps(ip - 1) - sqrt((dataset_ps(ip, 1) - dataset_ps(ip - 1, 1))^2 + (dataset_ps(ip, 2) - dataset_ps(ip - 1, 2))^2);
    s_ps(ip,2) = dataset_ps(ip,3);
    s_ps(ip,3) = dataset_ps(ip,4);
end
 
szn = [flip(s_ps(1:end,:)); s_ss];
end