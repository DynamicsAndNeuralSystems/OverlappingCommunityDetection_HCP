% Load the connectome
load ("RH.mat");

% Choose how many times to run Louvain
M = 100;
% Store community structure from each run
comms = zeros([180, M]);

for i=1:M
    comms(:,i) = community_louvain(RH);
end


% Compute pair-wise Jaccard coefficient for the above runs
JC = PairwiseJaccard(comms);

% Overlapping nodes acc to OSLOM
load('oslomcomm.mat');
overlapping = find(sum(oslomcomm>0, 2)>1);

% Average JC for overlapping nodes
mean(JC(overlapping))
% Average JC for all nodes
mean(JC)
