%-------------------------

% USER TO CHANGE: github repo base dir
github_repo_base_dir="/Users/abry4213/github/"

% Load stored benchmark network results
OCDA_HPC_repo_dir = fullfile(github_repo_base_dir, 'OverlappingCommunityDetection_HCP');
OCDA_comp_repo_dir = fullfile(github_repo_base_dir, 'OverlappingCommunityDetection');

% Add paths
addpath(genpath(OCDA_HPC_repo_dir));
addpath(genpath(OCDA_comp_repo_dir));

data_path = '/Users/abry4213/data/OCDA/';
RH = load(sprintf('%s/data/RH.mat', OCDA_HPC_repo_dir)).RH;

% Drop column/row 92, as it has all degree zero
RH(92, :) = [];  % Remove row 92
RH(:, 92) = [];  % Remove column 92


%% First-pass Louvain
% Run for 100 repeats
num_repeats = 100;
num_nodes = size(RH, 1);
gamma_val = 1.4;

% Sweep over gamma values from 1 to 1.5
gamma_list = 0.5:0.05:1.5;

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

%% 
% Looking specifically at gamma=1.05, let's compute the pairwise ENMI
% between each of the 100 seeds
target_gamma=1.05;
Louvain_at_target_gamma_res = readmatrix(sprintf('%s/Louvain_results/Louvain_Rubinov_assignments_100reps_gamma%0.2f.csv', data_path, target_gamma))


Louvain_at_target_gamma_NMI_res = zeros(10000, 3);
row_counter=1;

% Sweep over seeds for initialization from 1 to 100
seed_list = 1:1:100;

%%% ITERATE OVER TOLERANCES FOR EACH AS TARGET %%%
for target_seed_i = 1:length(seed_list)
    target_seed = seed_list(target_seed_i);
    target_seed_Louvain_res = Louvain_at_target_gamma_res(:,target_seed);

     %%% COMPARISON SEED RESULTS AT GIVEN TOLERANCE %%%
    for comparison_seed_i = 1:length(seed_list)
        comparison_seed = seed_list(comparison_seed_i);
        comparison_seed_Louvain_res = Louvain_at_target_gamma_res(:,comparison_seed);

        % Compute ENMI
        target_vs_comparison_seed_NMI = NMI_calc_vecs(target_seed_Louvain_res, comparison_seed_Louvain_res);

        % Add to table
        Louvain_at_target_gamma_NMI_res(row_counter, 1) = target_seed;
        Louvain_at_target_gamma_NMI_res(row_counter, 2) = comparison_seed;
        Louvain_at_target_gamma_NMI_res(row_counter, 3) = target_vs_comparison_seed_NMI;

        row_counter = row_counter+1;

    end
end

% Save the results
writematrix(Louvain_at_target_gamma_NMI_res, sprintf("%s/Louvain_results/NMI_results_sweep_across_seeds_%0.2f.csv", data_path, target_gamma));