% Takes as input the community assignments from multiple runs of a
% (non-overlapping) community detection algorithm on a given network. For
% each node, finds the average Jaccard coefficient of the communities it
% was assigned to over all pairs of runs. This serves as a metric of how
% robust the co-classifications of a given node are over multiple runs of
% community detection, and may help to identify nodes which don't fit
% cleanly into a single community ('overlapping' nodes).
%
% Input: 'CA' should be an NxM matrix, where N is the number of
% nodes in the network, and M is the number of times the community
% detection was run. Entry (i,j) of CA should be the community number to
% which node i was assigned on run j.
%
% Output: 'JC' is an Nx1 vector, with entry i giving the average Jaccard 
% coefficient over all run pairs for node i.
%
% Author: Sumeet Agarwal (sumeet@iitd.ac.in)
% Version: 26.09.2020

function JC = PairwiseJaccard(CA)
    N = size(CA,1); % Number of nodes
    M = size(CA,2); % Number of runs
    JC = zeros(N,1); % Vector of size N to store the average Jaccard coefficient values for all nodes
    for i=1:N
        for j=1:M
            for k=j+1:M
                cj = find(CA(:,j)==CA(i,j)); % Community of node i on run j
                ck = find(CA(:,k)==CA(i,k)); % Community of node i on run k
                jacc = length(intersect(cj,ck))/length(union(cj,ck)); % Jaccard coefficient of the two communities
                JC(i)= JC(i) + jacc/nchoosek(M,2); % Add on contribution to average Jaccard coefficient over all run pairs
            end
        end
    end    
end