%-------------------------
% ------------------------------------------------------
% Load benchmark network data
rng(0,'twister')

% USER TO CHANGE: github repo base dir
github_repo_base_dir="/Users/abry4213/github/"

% Load stored benchmark network results
OCDA_HPC_repo_dir = fullfile(github_repo_base_dir, 'OverlappingCommunityDetection_HCP');
OCDA_comp_repo_dir = fullfile(github_repo_base_dir, 'OverlappingCommunityDetection');

% storedBenchmark = fullfile(GiveMeFile('OCDA_toolbox'),'Computation','Conversions','networks');
networkDataFile = fullfile(OCDA_HPC_repo_dir, 'data', 'networks', 'network56.dat');
nodeLabelDataFile = fullfile(OCDA_HPC_repo_dir, 'data', 'networks', 'community56.dat');

% Run OCDAs on this benchmark
load(networkDataFile,'network56')

% Navigate to directory where OCDA algorithms are stored
cd(OCDA_comp_repo_dir);

%%
% Add path to OSLOM installation ()
addpath(fullfile(OCDA_comp_repo_dir, 'SourceCode', 'OSLOM'));
sourceCodePath = fullfile(OCDA_comp_repo_dir, 'SourceCode');
outputPath = fullfile(OCDA_HPC_repo_dir, 'data');
OLCD_Compute(network56, {'OSLOM', 'Shen', 'NNMF', 'Infomap', 'Jerry'}, true, nodeLabelDataFile, sourceCodePath, outputPath);

% rename output file
system(sprintf('mv %s/Computation_Result.mat %s/All_OCDA_Result_Representative_Benchmark.mat', outputPath, outputPath));

%%
% Visualize performance./../..
% cd ../
% cd Visualization_Module
OLCD_Visualization({'Benchmark', 'OSLOM', 'Shen', 'NNMF', 'Infomap','Jerry'}, sprintf('%s/All_OCDA_Result_Representative_Benchmark.mat', outputPath))


% ORIGINAL ADITI:
% %load benchmark network data
% generatesyntheticnetwork(1)
% network = load('weighted_networks/networks/network1.dat');

% cd /Users/aditijha/Desktop/communityDetection/
% %Run OCDAs on benchmark
% Computation(network, {'OSLOM', 'Clique', 'NNMF', 'Infomap','SLPA'}, 1, 'weighted_networks/communities/community1');

% % Visualize performance
% % cd ../
% % cd Visualization_Module
% Visualization({'Benchmark', 'OSLOM', 'Shen', 'NNMF', 'Infomap', 'SLPA'})
