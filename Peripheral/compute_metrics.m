% Code to generate 1000 synthetic networks [Doesn't have to be inside this]
num = 1000;
generatesyntheticnetwork(num);

% Matrix to store ENMI values
numMethods = 20;
ENMI_mat = zeros(num,numMethods);
F1_mat = zeros(num,numMethods);
Omega_mat = zeros(num,numMethods);


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
    % cd Visualization
    for j=1:size(Methods,1)
        % Compute ENMI
        ENMI_mat(i,j) = ENMI_calc(Final.Benchmark.Result, Final.(Methods{j}).Result);
        % Compute F1 score for overlapping nodes
        F1_mat(i,j) = F1_overlapcalc(Final.Benchmark.Result, Final.(Methods{j}).Result);
        % Compute Omega score
        % Omega_mat(i,j) = Omegaindex_calc(Final.Benchmark.Result, Final.(Methods{j}).Result);
    end
    % cd ..
    % cd Figures
end

save('ENMI_mat.mat','ENMI_mat');
save('F1_mat.mat','F1_mat');
% save('Omegaindex_mat.mat','OmegaIndex_mat');
