% This script makes group connectomes from HCP data

parcellation = 'HCPMMP1'; %'HCPMMP1' , 'custom200', custom500, aparcaseg
tract = 'iFOD2'; % FACT vs iFOD2
sift = 'SIFT2'; 
weight = 'standard'; % specify weight to be used:
% FA
% standard
% density
% MD
groupMatrixType = 'variance'; % specify type of group matrix
% 'variance' - Robert's method
% 'consistency' - VDH method
% 'lengthCV' - Misic's method
threshold = 0.90; % consistency threshold 
subjRem = 299; % subject in a cell 299 consistently has weird connectomes, so remove it

Conn = load(sprintf('%sANDfslatlas20_acpc_%s_%s_%s_structnets.mat', parcellation, tract, sift, weight));
Length = load(sprintf('%sANDfslatlas20_acpc_%s_%s_length_structnets.mat', parcellation, tract, sift));

coordinates = Conn.COG;
connectomes = Conn.ADJS;
distances = Length.ADJS;

% remove subjects with weird connectomes
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
    G = giveMeGroupAdj_variance(connectomes, 0.15); 
    % thisline can be changed to G = giveMeGroupAdj_variance(connectomes, dens); 
    % dens - will specify
    % the desired density of the connectome; by default the sensity will be
    % set to the average density of individual connectomes. 
elseif strcmp(groupMatrixType, 'consistency')
    [G,D] = giveMeGroupAdj_consistency(connectomes, distances, threshold);
end

