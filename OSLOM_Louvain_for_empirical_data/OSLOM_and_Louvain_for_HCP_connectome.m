% Define base repository
github_repo_base_dir=fileparts(pwd);

% Add all folders to path
addpath(genpath(github_repo_base_dir));

% Define connectome data path
connectome_data_path = fullfile(github_repo_base_dir, 'data', 'HCP_Connectome');

% Load in right-hemisphere structural connectome in HCP-MMP1 atlas
RH = load(fullfile(connectome_data_path, 'RH.mat')).RH;

% Drop column/row 92, as it has all degree zero
RH(92, :) = [];  % Remove row 92
RH(:, 92) = [];  % Remove column 92

%% Louvain

% Run Louvain with gamma=1 and seed=98
gamma=1.0;
Louvain_seed=98;

% Find number of nodes (should be 179 after dropping node 92)
num_nodes = size(RH, 1);

% Run Louvain
rng(Louvain_seed,'twister');
Louvain_community_affils = community_louvain(RH, gamma);

writematrix(Louvain_community_affils, ...
                sprintf('%s/data/louvain_community_assignments_BCT_gamma1_seed98.csv', github_repo_base_dir), ...
                Delimiter=',');

%% OSLOM
RH_data_copy='OSLOM_RH_copy.txt';

% Save the undirected adjacency matrix to a temporary txt file to be used
% for OSLOM
fid = fopen(RH_data_copy, 'w');
undirected_adj_mat = Mat2Undir(RH);
fprintf(fid,'%g\t%g\t%f\n', undirected_adj_mat');
fclose(fid);

tol=0.3;
OSLOM_seed=61;
numIters = 100;

% Define OSLOM code path
OSLOM_code_path=sprintf('%s/OCDA_evaluation/SourceCode/OSLOM', github_repo_base_dir);

% Define output file
OSLOM_output_file=sprintf('%s/data/OSLOM30_communities.csv', github_repo_base_dir);

% Run the command to start the OSLOM algorithm:
system(sprintf('%s/oslom_undir -f %s -w -r %g -t %f -seed %g', ...
    OSLOM_code_path, RH_data_copy, numIters, tol, OSLOM_seed));

% Remove unneeded files
system('rm -r OSLOM_RH_copy.txt_oslo_files');

fid = fopen('tp'); % Opens file
% Scans the text, creating a list, with 0s where it is the next module
A = textscan(fid, '%d', 'CommentStyle', '#', 'Delimiter', {' ', '\n'});
A = A{1}; % Sets it as the matrix, gets rid of the cell

% Closes and deletes the textfile
fclose(fid);
system('rm tp'); % Deletes the files
system(sprintf('rm %s', RH_data_copy));

numComms = sum(A==0);

% Creates the matrix for each tolerance:
commMat = zeros(num_nodes, numComms);

% Counter for communities
counter = 1;

for i = 1:length(A) % For the entire A list
    if ~A(i) == 0 % If it isn't 0
        commMat(A(i), counter) = 1; % Places a 1 in the community matrix
    else
        counter = counter + 1; % Else increases the counter
    end
end

for i = 1:num_nodes
    % Equally distributes the communities
    commMat(i, :) = commMat(i, :)/sum(commMat(i, :));
end

writematrix(commMat, OSLOM_output_file, Delimiter=',');