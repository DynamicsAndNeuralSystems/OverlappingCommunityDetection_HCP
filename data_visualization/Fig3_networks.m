%-------------------------
% ------------------------------------------------------
% Load benchmark network data
rng(0,'twister')

% Save the parent OCDA repo to add to path
parent_OCDA_repo = fileparts(pwd); % Get one level up

% Load stored benchmark network results
addpath(genpath(parent_OCDA_repo));

% Use the pre-computed network 122 as an example to visualize
networkDataFile = fullfile(parent_OCDA_repo, 'data', 'networks', 'network122.dat');
nodeLabelDataFile = fullfile(parent_OCDA_repo, 'data', 'communities', 'community122.dat');

% Run OCDAs on this benchmark
load(networkDataFile,'network122');


%%
% Navigate to directory where OCDA algorithms are stored
cd(parent_OCDA_repo);
output_benchmark_path = fullfile(parent_OCDA_repo, 'data');

% Add OSLOM to path specifically
addpath(fullfile(parent_OCDA_repo, 'OCDA_evaluation', 'SourceCode', 'OSLOM'));

% Define OCDAs
methods_list = {'Infomap', 'SLPA', 'OSLOM', 'Clique','NNMF'};

if ~isfile(sprintf('%s/All_OCDA_Result_Benchmark_Network122.mat',output_benchmark_path))
    sourceCodePath = fullfile(parent_OCDA_repo, 'OCDA_evaluation', 'SourceCode');
    OLCD_Compute(network122, methods_list, true, nodeLabelDataFile, sourceCodePath, output_benchmark_path);
    
    % rename output file
    system(sprintf('mv %s/Computation_Result.mat %s/All_OCDA_Result_Benchmark_Network122.mat', output_benchmark_path, output_benchmark_path));
else 
    load(sprintf('%s/All_OCDA_Result_Benchmark_Network122.mat',output_benchmark_path))
end

%%
% Plot the synthetic network adjacency matrix (sorted),
% ground-truth (simulated) overlapping community module assignment,
% along with community assignment across all algorithms per node
OLCD_Visualization(methods_list, ...
    sprintf('%s/All_OCDA_Result_Benchmark_Network122.mat', output_benchmark_path), ...
    "parula")
% print(gcf,'-vector','-dsvg','../plots/Exemplarbenchmark_heatmap.svg') % svg

%%
Methods = methods_list;
methodNames = cell(0); % Method names used
resultsFile = sprintf('%s/All_OCDA_Result_Benchmark_Network122.mat', output_benchmark_path);
cmap = 'parula';

% Load data
load(resultsFile);
fprintf('Loading data created on %s!\n', Final.Date);
numNodes = size(Final.Network, 1);

% Get the benchmark communities
benchComms = Final.Benchmark.Result;

% Benchmark sorting of nodes
I = Node_Sorter(benchComms); % Sorts the nodes into communities
disp('Nodes have been sorted into communities');