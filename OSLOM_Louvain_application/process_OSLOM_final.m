% Get working directory
repo_dir = fileparts(pwd);

% Define final tolerance, t_final
t_final = '0.30';

% Define final seed, 61
seed_final = 61;

% Define number of iterations
num_iters = 100;

% Define results path
oslomResPath = sprintf('%s/data/OSLOM_validation/', repo_dir);

% Defines the fileName for program to read
fileName = sprintf('%s/OSLOM_seed_%d_iters_%d_tol_%s.txt', oslomResPath, seed_final, num_iters, t_final);

fid = fopen(fileName); % Opens file
% Scans the text, creating a list, with 0s where it is the next module
A = textscan(fid, '%d', 'CommentStyle', '#', 'Delimiter', {' ', '\n'});
A = A{1}; % Sets it as the matrix, gets rid of the cell

% Closes the textfile
fclose(fid);

% There are 180 nodes in the Glasser right-hemisphere cortex
num_nodes=180;

% Find the number of communities
num_comms = sum(A==0);

% Creates the matrix for each tolerance:
OSLOM_comm_mat_wide = zeros(num_nodes, num_comms);

% Counter for communities
counter = 1;

for i = 1:length(A) % For the entire A list
    if ~A(i) == 0 % If it isn't 0
        OSLOM_comm_mat_wide(A(i), counter) = 1; % Places a 1 in the community matrix
    else
        counter = counter + 1; % Else increases the counter
    end
end

for i = 1:num_nodes
    % Equally distributes the communities
    OSLOM_comm_mat_wide(i, :) = OSLOM_comm_mat_wide(i, :)/sum(OSLOM_comm_mat_wide(i, :));
end

% Save OSLOM_comm_mat_wide
writematrix(OSLOM_comm_mat_wide, sprintf('%s/data/OSLOM30_final_module_assignments_wide_binary.csv', repo_dir));