% USER TO CHANGE: github repo base dir
github_repo_base_dir="/Users/abry4213/github/";

% Load stored benchmark network results
OCDA_HPC_repo_dir = fullfile(github_repo_base_dir, 'OverlappingCommunityDetection_HCP');
OCDA_comp_repo_dir = fullfile(github_repo_base_dir, 'OverlappingCommunityDetection');

% Add paths
addpath(OCDA_HPC_repo_dir);
addpath(OCDA_comp_repo_dir);

%% Load connectivity data and OSLOM results

% Right-hemisphere structural connectivity data -- 180x180
RH = load(sprintf('%s/data/RH.mat', OCDA_HPC_repo_dir)).RH;

% Drop index 92
RH(92, :) = [];  % Remove row 92
RH(:, 92) = [];  % Remove column 92

% OSLOM community assignments -- 180x7
OSLOM_final_results_wide = readmatrix(sprintf('%s/data/OSLOM30_final_module_assignments_wide_binary.csv', OCDA_HPC_repo_dir));

% Drop index 92
OSLOM_final_results_wide(92, :) = [];  % Remove row 92

% Load long version of OSLOM communities
OSLOM_final_results_long = readtable(sprintf('%s/data/OSLOM30_final_module_assignments.csv', OCDA_HPC_repo_dir), 'Format','%d %s %s %s');
original_community_labels = OSLOM_final_results_long.module;

% Find indices of overlapping nodes
overlapping_node_indices = find(sum(OSLOM_final_results_wide > 0, 2)>1); 

% Find the number of communities for each overlapping node
num_communities_per_node = sum(OSLOM_final_results_wide>0, 2);
num_communities_in_overlaps = num_communities_per_node(overlapping_node_indices);

%%
% Construct duplicate nodes
connectivity_mat = RH;
indices_to_repeat = overlapping_node_indices;
num_of_repeats = num_communities_in_overlaps - 1;
node_indices_original = OSLOM_final_results_long.node;

[RH_with_repeats, node_indices_with_repeats] = annie_overlapping_duplicates(connectivity_mat, indices_to_repeat, num_of_repeats, node_indices_original);

% Save results
save(sprintf('%s/data/RH_with_OSLOM_duplicates.mat', OCDA_HPC_repo_dir), "RH_with_repeats", "node_indices_with_repeats");



