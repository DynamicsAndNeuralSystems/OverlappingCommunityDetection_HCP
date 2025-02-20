% This script makes group connectomes from HCP data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define data path and imaging parameters
diffusion_data_path='../data/HCP_Connectome/';

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

% make a selected version of the group matrix
G = giveMeGroupAdj_variance(connectomes, dens); 

% Filter to the right hemisphere cortex
RH_full = G(hemiid==2, hemiid==2);
RH_cortex = RH_full(1:180, 1:180);

% Log-transform
RH_cortex_log = log(RH_cortex);