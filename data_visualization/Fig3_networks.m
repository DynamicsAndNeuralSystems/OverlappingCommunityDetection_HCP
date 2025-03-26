%-------------------------
% ------------------------------------------------------
% Load benchmark network data
rng(0,'twister')

% USER TO CHANGE: github repo base dir
github_repo_base_dir="/Users/abry4213/github/";

% Load stored benchmark network results
OCDA_HPC_repo_dir = fullfile(github_repo_base_dir, 'OverlappingCommunityDetection_HCP');
OCDA_comp_repo_dir = fullfile(github_repo_base_dir, 'OverlappingCommunityDetection');
addpath(genpath(OCDA_HPC_repo_dir));
addpath(genpath(OCDA_comp_repo_dir));

networkDataFile = fullfile(OCDA_HPC_repo_dir, 'data', 'networks', 'network56.dat');
nodeLabelDataFile = fullfile(OCDA_HPC_repo_dir, 'data', 'networks', 'community56.dat');

% Run OCDAs on this benchmark
load(networkDataFile,'network56')


%%
% Navigate to directory where OCDA algorithms are stored
cd(OCDA_comp_repo_dir);

outputPath = fullfile(OCDA_HPC_repo_dir, 'data');
addpath(fullfile(OCDA_comp_repo_dir, 'SourceCode', 'OSLOM'));

if ~isfile(sprintf('%s/All_OCDA_Result_Representative_Benchmark.mat',outputPath))
    % Add path to OSLOM installation ()
    addpath(fullfile(OCDA_comp_repo_dir, 'SourceCode', 'OSLOM'));
    sourceCodePath = fullfile(OCDA_comp_repo_dir, 'SourceCode');
    outputPath = fullfile(OCDA_HPC_repo_dir, 'data');
    OLCD_Compute(network56, {'OSLOM', 'Shen', 'NNMF', 'Infomap', 'Jerry'}, true, nodeLabelDataFile, sourceCodePath, outputPath);
    
    % rename output file
    system(sprintf('mv %s/Computation_Result.mat %s/All_OCDA_Result_Representative_Benchmark.mat', outputPath, outputPath));
else 
    load(sprintf('%s/All_OCDA_Result_Representative_Benchmark.mat',outputPath))
end


OLCD_Visualization({'Benchmark', 'OSLOM', 'Shen', 'NNMF', 'Infomap','Jerry'}, ...
    sprintf('%s/All_OCDA_Result_Representative_Benchmark.mat', outputPath), ...
    "parula")
% print(gcf,'-vector','-dsvg','../plots/Exemplarbenchmark_heatmap.svg') % svg