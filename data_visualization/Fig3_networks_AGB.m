%-------------------------
% ------------------------------------------------------
% Load benchmark network data
rng(0,'twister')

% Save the parent OCDA repo to add to path
parent_OCDA_repo = fileparts(pwd); % Get one level up

% Load stored benchmark network results
addpath(genpath(parent_OCDA_repo));

% Use the pre-computed network 803 as an example to visualize
networkDataFile = fullfile(parent_OCDA_repo, 'data', 'networks', 'network803.dat');
nodeLabelDataFile = fullfile(parent_OCDA_repo, 'data', 'communities', 'community803.dat');

% Run OCDAs on this benchmark
load(networkDataFile,'network803');


%%
% Navigate to directory where OCDA algorithms are stored
cd(parent_OCDA_repo);
output_benchmark_path = fullfile(parent_OCDA_repo, 'data', 'representative_benchmark');

% Add OSLOM to path specifically
addpath(fullfile(parent_OCDA_repo, 'OCDA_evaluation', 'SourceCode', 'OSLOM'));

% Define OCDAs
OCDA_methods_list = {'OSLOM', 'Clique', 'NNMF', 'Infomap', 'SLPA'};
methods_list = {'Benchmark', 'OSLOM', 'Clique', 'NNMF', 'Infomap', 'SLPA'};

if ~isfile(sprintf('%s/All_OCDA_Result_Benchmark_Network803.mat',output_benchmark_path))
    sourceCodePath = fullfile(parent_OCDA_repo, 'OCDA_evaluation', 'SourceCode');
    OLCD_Compute(network803, OCDA_methods_list, true, nodeLabelDataFile, sourceCodePath, output_benchmark_path);
    
    % rename output file
    system(sprintf('mv %s/Computation_Result.mat %s/All_OCDA_Result_Benchmark_Network803.mat', output_benchmark_path, output_benchmark_path));
else 
    load(sprintf('%s/All_OCDA_Result_Benchmark_Network803.mat',output_benchmark_path))
end
%%
% Save ground-truth community assignments
representative_benchmark_ground_truth_assignments = Final.Benchmark.Result;
writematrix(representative_benchmark_ground_truth_assignments, ...
    sprintf('%s/Ground_Truth_Comms_Network803.csv', output_benchmark_path), ...
        'Delimiter', ",")

% Extracting all method names
fields = fieldnames(Final);
fields = fields(4:end);
summary_stat_methods = fields;

for j=1:size(summary_stat_methods,1)
    this_OCDA = summary_stat_methods{j};
    % Load the community predictions with the given method
    this_OCDA_comm_res = Final.(summary_stat_methods{j}).Result;
    
    % Save the resulting OCDA community predictions matrix to a CSV
    writematrix(this_OCDA_comm_res, ...
        sprintf('%s/%s_Result_Benchmark_Network803.csv', output_benchmark_path, this_OCDA), ...
        'Delimiter', ",")
end

%%
% Plot the synthetic network adjacency matrix (sorted),
% ground-truth (simulated) overlapping community module assignment,
% along with community assignment across all algorithms per node
OLCD_Visualization(methods_list, ...
    sprintf('%s/All_OCDA_Result_Benchmark_Network803.mat', output_benchmark_path), ...
    "turbo")
print(gcf,'-vector','-dsvg','plots/benchmark_evaluation/Exemplarbenchmark_heatmap_network803.svg') % svg
