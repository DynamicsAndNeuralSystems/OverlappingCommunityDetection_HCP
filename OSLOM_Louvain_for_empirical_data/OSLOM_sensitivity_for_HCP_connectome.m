% Define base repository
github_repo_base_dir=fileparts(pwd);

% Add all folders to path
addpath(genpath(github_repo_base_dir));

% Define connectome data path
connectome_data_path = "/Users/abry4213/data/OCDA/HCP_Connectome/";

% Load in right-hemisphere structural connectome in HCP-MMP1 atlas
RH_cortex = load(fullfile(connectome_data_path, 'RH_cortex_0.15_CV_density.txt'));
LH_cortex = load(fullfile(connectome_data_path, 'LH_cortex_0.15_CV_density.txt'));

% Drop column/row 92, as it has all degree zero
RH_cortex(92, :) = [];  % Remove row 92
RH_cortex(:, 92) = [];  % Remove column 92

% Full cortex left/right
full_cortex = load(fullfile(connectome_data_path, 'Full_cortex_0.15_CV_density.txt'));

% Drop column/row 282, as it has all degree zero
full_cortex(282, :) = [];  % Remove row 92
full_cortex(:, 282) = [];  % Remove column 92

%% Right hemisphere only 
RH_data_copy='OSLOM_RH_copy.txt';
num_nodes = size(RH_cortex, 1);

% Save the undirected adjacency matrix to a temporary txt file to be used
% for OSLOM
fid = fopen(RH_data_copy, 'w');
undirected_adj_mat = Mat2Undir(RH_cortex);
fprintf(fid,'%g\t%g\t%f\n', undirected_adj_mat');
fclose(fid);

tol=0.3;
OSLOM_seed=61;
numIters = 100;

% Define OSLOM code path
OSLOM_code_path=sprintf('%s/OCDA_evaluation/SourceCode/OSLOM', github_repo_base_dir);

% Define output file
OSLOM_output_file=sprintf('%s/data/OSLOM30_robustness_RH_cortex_communities.csv', github_repo_base_dir);

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

%% Left hemisphere only 
LH_data_copy='OSLOM_LH_copy.txt';
num_nodes = size(LH_cortex, 1);

% Save the undirected adjacency matrix to a temporary txt file to be used
% for OSLOM
fid = fopen(LH_data_copy, 'w');
undirected_adj_mat = Mat2Undir(LH_cortex);
fprintf(fid,'%g\t%g\t%f\n', undirected_adj_mat');
fclose(fid);

tol=0.3;
OSLOM_seed=61;
numIters = 100;

% Define OSLOM code path
OSLOM_code_path=sprintf('%s/OCDA_evaluation/SourceCode/OSLOM', github_repo_base_dir);

% Define output file
OSLOM_output_file=sprintf('%s/data/OSLOM30_robustness_LH_cortex_communities.csv', github_repo_base_dir);

% Run the command to start the OSLOM algorithm:
system(sprintf('%s/oslom_undir -f %s -w -r %g -t %f -seed %g', ...
    OSLOM_code_path, LH_data_copy, numIters, tol, OSLOM_seed));

% Remove unneeded files
system('rm -r OSLOM_LH_copy.txt_oslo_files');

fid = fopen('tp'); % Opens file
% Scans the text, creating a list, with 0s where it is the next module
A = textscan(fid, '%d', 'CommentStyle', '#', 'Delimiter', {' ', '\n'});
A = A{1}; % Sets it as the matrix, gets rid of the cell

% Closes and deletes the textfile
fclose(fid);
system('rm tp'); % Deletes the files
system(sprintf('rm %s', LH_data_copy));

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

%% Both hemispheres together
full_data_copy='OSLOM_full_copy.txt';
num_nodes = size(full_cortex, 1);

% Save the undirected adjacency matrix to a temporary txt file to be used
% for OSLOM
fid = fopen(full_data_copy, 'w');
undirected_adj_mat = Mat2Undir(full_cortex);
fprintf(fid,'%g\t%g\t%f\n', undirected_adj_mat');
fclose(fid);

tol=0.3;
OSLOM_seed=61;
numIters = 100;

% Define OSLOM code path
OSLOM_code_path=sprintf('%s/OCDA_evaluation/SourceCode/OSLOM', github_repo_base_dir);

% Define output file
OSLOM_output_file=sprintf('%s/data/OSLOM30_robustness_full_cortex_communities.csv', github_repo_base_dir);

% Run the command to start the OSLOM algorithm:
system(sprintf('%s/oslom_undir -f %s -w -r %g -t %f -seed %g', ...
    OSLOM_code_path, full_data_copy, numIters, tol, OSLOM_seed));

% Remove unneeded files
system('rm -r OSLOM_full_copy.txt_oslo_files');

fid = fopen('tp'); % Opens file
% Scans the text, creating a list, with 0s where it is the next module
A = textscan(fid, '%d', 'CommentStyle', '#', 'Delimiter', {' ', '\n'});
A = A{1}; % Sets it as the matrix, gets rid of the cell

% Closes and deletes the textfile
fclose(fid);
system('rm tp'); % Deletes the files
system(sprintf('rm %s', full_data_copy));

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

