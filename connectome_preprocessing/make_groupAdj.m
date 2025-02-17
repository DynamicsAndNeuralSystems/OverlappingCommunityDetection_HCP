% This script makes group connectomes from HCP data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define data path and imaging parameters
diffusion_data_path='../data/HCP_Connectome/';

% Parcellation: options are 'HCPMMP1' , 'custom200', 'custom500',
% 'aparc+aseg'
% We use 'HCPMMP1' to work with the Glasser (2016) 360-region atlas
parcellation = 'HCPMMP1';

% Tractography type: FACT or iFOD2; we use iFOD2 here
tract = 'iFOD2';

% SIFT type: SIFT or SIFT2; we use SIFT2 here
sift = 'SIFT2'; 

% Weighting to be used: FA, standard, density, or MD; we use standard
weight = 'standard';

% How to generate group adjacency matrix:
% 'variance' - Robert's method
% 'consistency' - VDH method
% 'lengthCV' - Misic's method
% groupMatrixType = 'variance';
groupMatrixType = 'variance';
threshold = 0.9;
dens = 0.15;

% Load data
SC_data = load(sprintf('%s/%sANDfslatlas20_acpc_%s_%s_%s_structnets.mat', diffusion_data_path, parcellation, tract, sift, weight));
length_data = load(sprintf('%s/%sANDfslatlas20_acpc_%s_%s_length_structnets.mat', diffusion_data_path, parcellation, tract, sift));

% Extract coordinates, connectomes, and distances
coordinates = SC_data.COG;
connectomes = SC_data.ADJS;
distances = length_data.ADJS;

% remove subjects with weird connectomes
subjRem = 299; % subject in a cell 299 consistently has weird connectomes, so remove it
connectomes(subjRem) = []; 
coordinates(subjRem) = []; 
distances(subjRem) = []; 

% if there are any emtpy cells, remove them as well
% normally there will be no empty cells, but in case some subjects failed
% connectome generation
connectomes = connectomes(~cellfun('isempty',connectomes));
coordinates = coordinates(~cellfun('isempty',coordinates));
distances = distances(~cellfun('isempty',distances));

% calculate average distances between ROIs as an average of all subjects
numNodes = size(coordinates{1},1);
numSubj = size(coordinates,2);

% make vectors for hemispheres for different parcellations based on the
% number of nodes (first half is always left, second half is always right)
hemiid = zeros(numNodes,1);
hemiid(1:numNodes/2) = 1;
hemiid(numNodes/2+1:numNodes) = 2;

dist = zeros(numNodes, numNodes, numSubj);
adjMatr = zeros(numNodes, numNodes, numSubj);

for s=1:numSubj
    dist(:,:,s) = pdist2(coordinates{s}, coordinates{s});
    adjMatr(:,:,s) = connectomes{s};
end

meanDist = mean(dist,3); 

% make a selected version of the group matrix
if strcmp(groupMatrixType, 'lengthCV')
    G = fcn_group_average(adjMatr,meanDist,hemiid);
elseif strcmp(groupMatrixType, 'variance')
    % G = giveMeGroupAdj_variance(connectomes, 0.15); 
    G = giveMeGroupAdj_variance(connectomes, dens); 
    % thisline can be changed to G = giveMeGroupAdj_variance(connectomes, dens); 
    % dens - will specify
    % the desired density of the connectome; by default the sensity will be
    % set to the average density of individual connectomes. 
elseif strcmp(groupMatrixType, 'consistency')
    [G,D] = giveMeGroupAdj_consistency(connectomes, distances, threshold);
end

% Filter to the right hemisphere cortex
RH_full = G(hemiid==2, hemiid==2);
RH_cortex = RH_full(1:180, 1:180);