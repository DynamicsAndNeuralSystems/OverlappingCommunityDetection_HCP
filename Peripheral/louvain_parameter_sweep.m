%-------------------------

% USER TO CHANGE: github repo base dir
github_repo_base_dir="/Users/abry4213/github/"

% Load stored benchmark network results
OCDA_HPC_repo_dir = fullfile(github_repo_base_dir, 'OverlappingCommunityDetection_HCP');
OCDA_comp_repo_dir = fullfile(github_repo_base_dir, 'OverlappingCommunityDetection');

% Add paths
addpath(OCDA_HPC_repo_dir);
addpath(OCDA_comp_repo_dir);

data_path = '/Users/abry4213/data/OCDA/';
RH = load(sprintf('%s/RH.mat', data_path)).RH;

%% First-pass Louvain
% Run for 100 repeats
num_repeats = 100;
num_nodes = size(RH, 1);
gamma_val = 1.4;

% Sweep over gamma values from 1 to 1.5
gamma_list = 1:0.05:1.5;

for gamma_list_i = 1:length(gamma_list)
    gamma_val = gamma_list(gamma_list_i);
    community_affils = zeros(num_nodes, num_repeats);
    for repeat_i = 1:num_repeats
        % Run Louvain
        rng(repeat_i,'twister');
        repeat_louvain_res = community_louvain(RH, gamma_val);
        community_affils(:,repeat_i) = repeat_louvain_res;
    end
    
    writematrix(community_affils, ...
                sprintf('%s/Louvain_results/Louvain_Rubinov_assignments_100reps_gamma%0.2f.csv', data_path, gamma_val), ...
                Delimiter=',');
end
