%% Load data 
data_path = '/Users/abry4213/github/OverlappingCommunityDetection_HCP/Results/';
RH = load(sprintf('%s/RH.mat', data_path)).RH;

%%
% Add OCD code to path
OCD_code_path = '/Users/abry4213/github/OverlappingCommunityDetection/';
addpath(OCD_code_path);

% Create undirected adjacency matrix from the RH connectome
undirected_adj_mat = Mat2Undir(RH);

% Find the number of nodes
num_nodes = size(RH, 1);

% Use 100 iterations
num_iters = 100;

% Sweep through tolerance values from 0.15 to 0.35
Tol = 0.14:0.01:0.35;

% Define OSLOM code path
OSLOM_code_path ='/Users/abry4213/Downloads/OSLOM2';

%% 
% Save the undirected adjacency matrix to a temporary txt file
fid = fopen('tempData.txt', 'w');
fprintf(fid,'%g\t%g\t%f\n', undirected_adj_mat');
fclose(fid);

%%
%-------------------------------------------------------------------------------
% Run OSLOM
%-------------------------------------------------------------------------------
num_tol = length(Tol);
fprintf(1,'Running OSLOM on the network across %u tolerances\n',num_tol);
for t = 1:num_tol
    tol = Tol(t);
    % Run the command to start the OSLOM algorithm:
    system(sprintf('%s/oslom_undir -f tempData.txt -w -r %g -t %f', OSLOM_code_path, num_iters, tol));

    % Move/rename the important data to base location:
    system(sprintf('mv tp %s/OSLOM_validation/OSLOM_tol_%g.txt', data_path, tol));
    
    % Remove unneeded files:
    system('rm -r tempData.txt_oslo_files');

end

%%
% Process the final results 

% Use 0.3 as the final tolerance
Tol_final = 0.3;

% Creates the cells for the output
Output = cell(1, length(Tol_final));

for tol = Tol_final
    % Defines the fileName for program to read
    fileName = sprintf('%s/OSLOM_validation/OSLOM_tol_%g.txt', data_path, tol);

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

    Output{Tol_final == tol} = commMat; % Sets it as the output
end

%%
oslomcomm_mine = Output{1};