% Code to generate 1000 synthetic networks [Doesn't have to be inside this]
% The code then computes extended normalized mutual information (ENMI) for
% each of the 20 methods on each of the 1000 synthetic networks
num = 1000;
generatesyntheticnetwork(num);

% Matrix to store ENMI values
numMethods = 22;
ENMI_mat = zeros(num,numMethods);

for i = 1:num
    % Load each synthetic networks previously generated
    benchpath = ['weighted_networks/networks/network', int2str(i), '.dat'];
    benchmark = load(benchpath);
    commpath = ['weighted_networks/communities/community', int2str(i), '.dat'];

    % Run OCDAs
    Final = Computation(benchmark, {'OSLOM', 'Infomap','SLPA','Clique','NNMF'}, 1, commpath);
    fields = fieldnames(Final);
    fields = fields(4:end);

    % Extracting all method names
    Methods = fields;
    
    % Compute ENMI for each method
    for j=1:size(Methods,1)
        % Compute ENMI
        ENMI_mat(i,j) = ENMI_calc(Final.Benchmark.Result, Final.(Methods{j}).Result);
    end
end

save('ENMI_mat.mat','ENMI_mat');
