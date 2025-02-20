% USER TO CHANGE: github repo base dir
github_repo_base_dir="/Users/abry4213/github/";

% Load stored benchmark network results
OCDA_HPC_repo_dir = fullfile(github_repo_base_dir, 'OverlappingCommunityDetection_HCP');
OCDA_comp_repo_dir = fullfile(github_repo_base_dir, 'OverlappingCommunityDetection');

% Add paths
addpath(OCDA_HPC_repo_dir);
addpath(OCDA_comp_repo_dir);

data_path = '/Users/abry4213/data/OCDA/';

%% Load connectivity data
RH = load(sprintf('%s/RH.mat', data_path)).RH;

% Drop index 92
RH(92, :) = [];  % Remove row 92
RH(:, 92) = [];  % Remove column 92

% Find indices of overlapping nodes
overlapping_node_indices = find(sum(OSLOM_final_results_wide > 0, 2)>1); 

%% Construct the overlapping duplicates
[RH_with_duplicates,commLabelsNew,iOverlappingNew,nodeListTracking] = constructOverlappingDuplicates(RH,overlapping_node_indices,OSLOM_final_results_wide);
save(sprintf('%s/data/RH_new.mat', OCDA_HPC_repo_dir), "iOverlappingNew", "RH_with_duplicates", "nodeListTracking", "commLabelsNew");
