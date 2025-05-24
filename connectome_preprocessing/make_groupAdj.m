% This script makes group connectomes from HCP data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define data path and imaging parameters
% diffusion_data_path='../data/HCP_Connectome/';
diffusion_data_path='/Users/abry4213/data/OCDA/HCP_Connectome/';

% Add Brain Connectivity Toolbox to path
addpath('/Users/abry4213/Documents/MATLAB/BCT/2019_03_03_BCT/');

% We use 'HCPMMP1' to work with the Glasser (2016) 360-region atlas
parcellation = 'HCPMMP1';

% SIFT type: SIFT2 
sift = 'SIFT2'; 

% Weighting to be used: standard
weight = 'standard';

% Tractography type: iFOD2
tract = 'iFOD2';

% How to generate group adjacency matrix:
% 'variance' - Robert's method
% density of 0.15 means keeping edges where the weight CV is in the bottom
% 15%
groupMatrixType = 'variance';
dens = 0.15;

% Load data
SC_data = load(sprintf('%s/%sANDfslatlas20_acpc_%s_%s_%s_structnets.mat', diffusion_data_path, parcellation, tract, sift, weight));

% Extract coordinates, connectomes, and distances
coordinates = SC_data.COG;
connectomes = SC_data.ADJS;

% remove subjects with weird connectomes
subjRem = 299; % subject in a cell 299 consistently has weird connectomes, so remove it
connectomes(subjRem) = []; 
coordinates(subjRem) = []; 
% if there are any emtpy cells, remove them as well
% normally there will be no empty cells, but in case some subjects failed
% connectome generation
connectomes = connectomes(~cellfun('isempty',connectomes));
coordinates = coordinates(~cellfun('isempty',coordinates));

% make vectors for hemispheres for different parcellations based on the
% number of nodes (first half is always left, second half is always right)
numNodes = size(coordinates{1},1);
hemiid = zeros(numNodes,1);
hemiid(1:numNodes/2) = 1;
hemiid(numNodes/2+1:numNodes) = 2;

%%
% No thresholding -- check (1) what the unthresholded connectome looks
% like, (2) what the CV is for all edges
[group_mat_no_thresh, group_mat_CV] = giveMeGroupAdj_variance(connectomes, 1); 

% Log-transform
group_mat_no_thresh = log(group_mat_no_thresh);
group_mat_no_thresh(isinf(group_mat_no_thresh)) = 0 

% Save the unthresholded matrix and CV
writematrix(group_mat_no_thresh, ...
    sprintf('%s/Full_brain_adj_mat_no_thresholding.txt', diffusion_data_path))
writematrix(group_mat_CV, ...
    sprintf('%s/Full_brain_edge_CV_no_thresholding.txt', diffusion_data_path))

% Also filter to just the right-hemisphere cortex for sensitivity analysis
RH_full_no_thresh = group_mat_no_thresh(hemiid==2, hemiid==2);
RH_cortex_no_thresh = RH_full_no_thresh(1:180, 1:180);

RH_full_CV_no_thresh = group_mat_CV(hemiid==2, hemiid==2);
RH_cortex_CV_no_thresh = RH_full_CV_no_thresh(1:180, 1:180);

% Save the unthresholded right cortex matrix and CV
writematrix(RH_cortex_no_thresh, ...
    sprintf('%s/RH_cortex_adj_mat_no_thresholding.txt', diffusion_data_path))
writematrix(RH_cortex_CV_no_thresh, ...
    sprintf('%s/RH_cortex_edge_CV_no_thresholding.txt', diffusion_data_path))


%% Sweep across different thresholds to compare resulting connectome
% Main analysis uses dens = 0.15 
for dens = [0.15, 0.3, 0.45, 0.6]
    [G, G_CV] = giveMeGroupAdj_variance(connectomes, dens); 

    % Log-transform 
    G = log(G);

    % Set inf to 0
    G(isinf(G)) = 0;     % find inf values and replace with 0

    % Filter to left/right cortex
    cortex_indices = [1:180, 191:370];
    G_cortex = G(cortex_indices,cortex_indices);

    % Save the whole-brain connectome at this density
    writematrix(G_cortex, sprintf('%s/Full_cortex_%.2f_CV_density.txt', diffusion_data_path, dens))

    % Subset down to left hemisphere
    LH_full = G(hemiid==1, hemiid==1);
    LH_cortex = LH_full(1:180, 1:180);
    writematrix(LH_cortex, sprintf('%s/LH_cortex_%.2f_CV_density.txt', diffusion_data_path, dens))

    % Subset down to right hemisphere
    RH_full = G(hemiid==2, hemiid==2);
    RH_cortex = RH_full(1:180, 1:180);
    writematrix(RH_cortex, sprintf('%s/RH_cortex_%.2f_CV_density.txt', diffusion_data_path, dens))

end
