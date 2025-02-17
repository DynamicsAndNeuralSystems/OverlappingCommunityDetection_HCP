%% Load data 
data_path = '/Users/abry4213/data/OCDA/';
RH = load(sprintf('%s/RH.mat', data_path)).RH;

% Add OCD code to path
OCD_code_path = '/Users/abry4213/github/OverlappingCommunityDetection/';
addpath(OCD_code_path);

%%
% There are 180 nodes in the network
num_nodes=180;

% Sweep over p-values from 0.01 to 1
tol_list = 0.01:0.01:1;

% Sweep over seeds for initialization from 1 to 100
seed_list = 1:1:100;

%%

% Create undirected adjacency matrix from the RH connectome
undirected_adj_mat = Mat2Undir(RH);
% Save the undirected adjacency matrix to a temporary txt file to be used
% for OSLOM
fid = fopen(sprintf('%s/RH.txt', data_path), 'w');
fprintf(fid,'%g\t%g\t%f\n', undirected_adj_mat');
fclose(fid);

%% Analyzing ENMI for each p-value across seeds
ENMI_results_by_tol_sweep_across_seeds = zeros(100000, 4);
row_counter=1;

for target_tol_i = 1:length(tol_list)
    target_tol = tol_list(target_tol_i);

    %%% ITERATE OVER TOLERANCES FOR EACH AS TARGET %%%
    for target_seed_i = 1:length(seed_list)
        target_seed = seed_list(target_seed_i);
        
        % Load in results for this target seed at this target tolerance
        OSLOM_target_tol_target_seed_file = sprintf('%s/OSLOM_results/OSLOM_seed_%d_iters_100_tol_%0.2f.txt', ...
                                    data_path, target_seed, target_tol);

        fid_target_tol_target_seed = fopen(OSLOM_target_tol_target_seed_file);
        % Scans the text, creating a list, with 0s where it is the next module
        A_target_tol_target_seed = textscan(fid_target_tol_target_seed, '%d', 'CommentStyle', '#', 'Delimiter', {' ', '\n'});
        A_target_tol_target_seed = A_target_tol_target_seed{1}; % Sets it as the matrix, gets rid of the cell

        % Close file
        fclose(fid_target_tol_target_seed);
        
        num_modules_target_tol_target_seed = sum(A_target_tol_target_seed==0);
        
        % Creates the matrix for each tolerance:
        OSLOM_target_tol_target_seed_res = zeros(num_nodes, num_modules_target_tol_target_seed);
        
        % Counter for communities
        counter = 1;
        
        for A_i = 1:length(A_target_tol_target_seed) % For the entire A list
            if ~A_target_tol_target_seed(A_i) == 0 % If it isn't 0
                OSLOM_target_tol_target_seed_res(A_target_tol_target_seed(A_i), counter) = 1; % Places a 1 in the community matrix
            else
                counter = counter + 1; % Else increases the counter
            end
        end
        
        for node_i = 1:num_nodes
            % Equally distributes the communities
            OSLOM_target_tol_target_seed_res(node_i, :) = OSLOM_target_tol_target_seed_res(node_i, :)/sum(OSLOM_target_tol_target_seed_res(node_i, :));
        end

        %%% COMPARISON SEED RESULTS AT GIVEN TOLERANCE %%%
        for comparison_seed_i = 1:length(seed_list)
            comparison_seed = seed_list(comparison_seed_i);
            % Load in results for this target seed at this target tolerance
            OSLOM_target_tol_comparison_seed_file = sprintf('%s/OSLOM_results/OSLOM_seed_%d_iters_100_tol_%0.2f.txt', ...
                                        data_path, comparison_seed, target_tol);
    
            fid_target_tol_comparison_seed = fopen(OSLOM_target_tol_comparison_seed_file);
            % Scans the text, creating a list, with 0s where it is the next module
            A_target_tol_comparison_seed = textscan(fid_target_tol_comparison_seed, '%d', 'CommentStyle', '#', 'Delimiter', {' ', '\n'});
            A_target_tol_comparison_seed = A_target_tol_comparison_seed{1}; % Sets it as the matrix, gets rid of the cell

            fclose(fid_target_tol_comparison_seed);
            
            num_modules_target_tol_comparison_seed = sum(A_target_tol_comparison_seed==0);
            
            % Creates the matrix for each tolerance:
            OSLOM_target_tol_comparison_seed_res = zeros(num_nodes, num_modules_target_tol_comparison_seed);
            
            % Counter for communities
            counter = 1;
            
            for A_i = 1:length(A_target_tol_comparison_seed) % For the entire A list
                if ~A_target_tol_comparison_seed(A_i) == 0 % If it isn't 0
                    OSLOM_target_tol_comparison_seed_res(A_target_tol_comparison_seed(A_i), counter) = 1; % Places a 1 in the community matrix
                else
                    counter = counter + 1; % Else increases the counter
                end
            end
            
            for node_i = 1:num_nodes
                % Equally distributes the communities
                OSLOM_target_tol_comparison_seed_res(node_i, :) = OSLOM_target_tol_comparison_seed_res(node_i, :)/sum(OSLOM_target_tol_comparison_seed_res(node_i, :));
            end

            % COMPUTE ENMI TO TARGET SEED
            ENMI_target_tol_target_seed = ENMI_calc(OSLOM_target_tol_target_seed_res, OSLOM_target_tol_comparison_seed_res);

            % Append to results
            ENMI_results_by_tol_sweep_across_seeds(row_counter, 1) = target_tol;
            ENMI_results_by_tol_sweep_across_seeds(row_counter, 2) = target_seed;
            ENMI_results_by_tol_sweep_across_seeds(row_counter, 3) = comparison_seed;
            ENMI_results_by_tol_sweep_across_seeds(row_counter, 4) = ENMI_target_tol_target_seed;

            % Increase row_counter
            row_counter = row_counter + 1;

        end

    end
end

% Save the results
writematrix(ENMI_results_by_tol_sweep_across_seeds, sprintf("%s/OSLOM_results/ENMI_results_by_tol_sweep_across_seeds.csv", data_path))

