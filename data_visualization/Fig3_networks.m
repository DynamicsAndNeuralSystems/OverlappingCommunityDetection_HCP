%-------------------------
% ------------------------------------------------------
% Load benchmark network data
rng(0,'twister')

% Save the parent OCDA repo to add to path
parent_OCDA_repo = fileparts(pwd); % Get one level up

% Load stored benchmark network results
addpath(genpath(parent_OCDA_repo));

% Use the pre-computed network 56 as an example to visualize
networkDataFile = fullfile(parent_OCDA_repo, 'data', 'networks', 'network56.dat');
nodeLabelDataFile = fullfile(parent_OCDA_repo, 'data', 'networks', 'community56.dat');

% Run OCDAs on this benchmark
load(networkDataFile,'network56');


%%
% Navigate to directory where OCDA algorithms are stored
cd(parent_OCDA_repo);
output_benchmark_path = fullfile(parent_OCDA_repo, 'data');

% Add OSLOM to path specifically
addpath(fullfile(parent_OCDA_repo, 'evaluate_OCDAs_on_benchmark', 'SourceCode', 'OSLOM'));

if ~isfile(sprintf('%s/All_OCDA_Result_Representative_Benchmark.mat',output_benchmark_path))
    sourceCodePath = fullfile(parent_OCDA_repo, 'evaluate_OCDAs_on_benchmark', 'SourceCode');
    OLCD_Compute(network56, {'OSLOM', 'Shen', 'NNMF', 'Infomap', 'Jerry'}, true, nodeLabelDataFile, sourceCodePath, output_benchmark_path);
    
    % rename output file
    system(sprintf('mv %s/Computation_Result.mat %s/All_OCDA_Result_Representative_Benchmark.mat', output_benchmark_path, output_benchmark_path));
else 
    load(sprintf('%s/All_OCDA_Result_Representative_Benchmark.mat',output_benchmark_path))
end

% Plot the synthetic network adjacency matrix (sorted),
% ground-truth (simulated) overlapping community module assignment,
% along with community assignment across all algorithms per node
OLCD_Visualization({'Benchmark', 'OSLOM', 'Shen', 'NNMF', 'Infomap','Jerry'}, ...
    sprintf('%s/All_OCDA_Result_Representative_Benchmark.mat', output_benchmark_path), ...
    "turbo")
% print(gcf,'-vector','-dsvg','../plots/Exemplarbenchmark_heatmap.svg') % svg