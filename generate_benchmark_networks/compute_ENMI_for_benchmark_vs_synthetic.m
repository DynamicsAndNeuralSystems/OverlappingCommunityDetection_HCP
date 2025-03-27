% Generate 1000 synthetic networks using generate_synthetic_networks.m,
% Then compute extended normalized mutual information (ENMI) for
% each of the 22 atatlmethods on each of the 1000 synthetic networks
% Two arguments to generate_synthetic_networks: number of networks to generate and the directory to save the networks
num_networks = 1000;
network_output_dir = [fileparts(pwd), '/data/'];
generate_synthetic_networks(num_networks, network_output_dir);

% Add computation code to path
OCDA_evaluation_path = sprintf('%s/OCDA_evaluation', fileparts(pwd));
addpath(genpath(fullfile(OCDA_evaluation_path, 'Computation')));
addpath(genpath(fullfile(OCDA_evaluation_path, 'SourceCode')));

sourceCodePath = fullfile(OCDA_evaluation_path, 'SourceCode');

% Matrix to store ENMI values
numMethods = 22;
ENMI_mat = zeros(num_networks,numMethods);

for i = 1:10
    % Load each synthetic networks previously generated
    benchpath = sprintf('%s/networks/network%s.dat', network_output_dir, int2str(i));
    benchmark = load(benchpath);
    commpath = sprintf('%s/communities/community%s.dat', network_output_dir, int2str(i));

    % Run OCDAs: networkAdj,methodList,isBenchmark,benchmarkFileName,sourceCodePath, outputPath
    methods_list = {'SLPA', 'Infomap', 'OSLOM', 'Clique','NNMF'};
    Final = OLCD_Compute(benchmark, methods_list, 1, ...
        commpath, sourceCodePath, sprintf('%s/ENMI_Results', network_output_dir));
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

save(sprintf('%s/ENMI_Results/ENMI_mat.mat', network_output_dir),'ENMI_mat');
