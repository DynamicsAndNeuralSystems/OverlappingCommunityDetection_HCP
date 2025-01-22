%% Load data 
data_path = '/Users/abry4213/github/OverlappingCommunityDetection_HCP/Results/';
RH = load(sprintf('%s/RH.mat', data_path)).RH;

% Add OCD code to path
OCD_code_path = '/Users/abry4213/github/OverlappingCommunityDetection/';
addpath(OCD_code_path);

% Create undirected adjacency matrix from the RH connectome
undirected_adj_mat = Mat2Undir(RH);

% Find the number of nodes
num_nodes = size(RH, 1);

% Use 100 iterations
num_iters = 100;

% Set the seed to 127
seed_num = 127;

% Sweep through tolerance values from 0.15 to 0.35
Tol = 0.01:0.01:1;

% Define OSLOM code path
OSLOM_code_path ='/Users/abry4213/Downloads/OSLOM2';

% Save the undirected adjacency matrix to a temporary txt file
fid = fopen('tempData.txt', 'w');
fprintf(fid,'%g\t%g\t%f\n', undirected_adj_mat');
fclose(fid);

%-------------------------------------------------------------------------------
% Run OSLOM
%-------------------------------------------------------------------------------

num_tol = length(Tol);
fprintf(1,'Running OSLOM on the network across %u tolerances\n',num_tol);
for t = 1:num_tol
    tol = Tol(t);
    % Check if output file doesn't yet exist
    tol_output_file = sprintf('%s/OSLOM_validation/OSLOM_seed_%d_tol_%g.txt', data_path, seed_num, tol);
    if ~isfile(tol_output_file)
        % Run the command to start the OSLOM algorithm:
        system(sprintf('%s/oslom_undir -f tempData.txt -w -r %g -t %f -seed %d', OSLOM_code_path, num_iters, tol, seed_num));
         
        % Move/rename the important data to base location:
        system(sprintf('mv tp %s', tol_output_file));
        
        % Remove unneeded files:
        system('rm -r tempData.txt_oslo_files');
    end
end

%% Find the number of communities as a result of each tolerance level
num_comms = zeros(length(Tol),2);

for i = 1:length(Tol)
    tol = Tol(i);

    % Defines the fileName for program to read
    fileName = sprintf('%s/OSLOM_validation/OSLOM_seed_%d_tol_%g.txt', data_path, seed_num, tol);
    fid = fopen(fileName); % Opens file

    % Scans the text, creating a list, with 0s where it is the next module
    A = textscan(fid, '%d', 'CommentStyle', '#', 'Delimiter', {' ', '\n'});
    A = A{1}; % Sets it as the matrix, gets rid of the cell

    % Closes and deletes the textfile
    fclose(fid);

    % Finds the total number of communities detected at current p-value
    % threshold, tol
    numComms = sum(A==0);

    % Append the number of detected communities to num_comms
    num_comms(i, 1) = tol;
    num_comms(i, 2) = numComms;

end

% save the num_comms matrix to a table
dlmwrite(sprintf('%s/OSLOM_validation/OSLOM_seed_%d_num_comms.txt', data_path, seed_num), num_comms, 'delimiter','\t')

%%
% Process the final results, using 0.3 as the final tolerance
Tol_final = 0.3;

% Creates the cells for the output
Output = cell(1, length(Tol_final));

% Defines the fileName for program to read
fileName = sprintf('%s/OSLOM_validation/OSLOM_seed_%d_tol_%g.txt', data_path, seed_num, Tol_final);

fid = fopen(fileName); % Opens file
% Scans the text, creating a list, with 0s where it is the next module
A = textscan(fid, '%d', 'CommentStyle', '#', 'Delimiter', {' ', '\n'});
A = A{1}; % Sets it as the matrix, gets rid of the cell

% Closes and deletes the textfile
fclose(fid);

% command = sprintf('rm %s',fileName);
% system(command); % Deletes the file

numComms = sum(A==0);

% Creates the matrix for each tolerance:
oslomcomm_AGB = zeros(num_nodes, numComms);

% Counter for communities
counter = 1;

for i = 1:length(A) % For the entire A list
    if ~A(i) == 0 % If it isn't 0
        oslomcomm_AGB(A(i), counter) = 1; % Places a 1 in the community matrix
    else
        counter = counter + 1; % Else increases the counter
    end
end

for i = 1:num_nodes
    % Equally distributes the communities
    oslomcomm_AGB(i, :) = oslomcomm_AGB(i, :)/sum(oslomcomm_AGB(i, :));
end

save(sprintf('%s/OSLOM_validation/oslomcomm_seed_%d_AGB.mat', data_path, seed_num), 'oslomcomm_AGB');
