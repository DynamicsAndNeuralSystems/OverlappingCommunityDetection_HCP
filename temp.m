%%
load("Results/RH.mat");
load("Results/RH_new.mat");
load("Results/louvaincomm_Rubinov.mat");
load("Results/oslomcomm.mat");
iOverlapping = find(sum(oslomcomm > 0, 2)==2); % list of overlapping nodes obtained from OSLOM

%%
adj=RH;
commLabels=oslomcomm;

numNodes = length(adj);
numOverlapping = length(iOverlapping);
newNodes = numNodes + numOverlapping;

% Duplicate the main adj
adjNew = zeros(newNodes,newNodes);
adjNew(1:numNodes,1:numNodes) = adj;

% Create list with node indices from original oslom_mat
nodeListTracking = zeros(newNodes,1);
nodeListTracking(1:numNodes) = 1:numNodes;
nodeListTracking(numNodes+1:newNodes) = iOverlapping;

% Add the overlapping nodes as additional nodes
rangeOrig = 1:numNodes;
rangeNew = numNodes+1:numNodes+numOverlapping;
adjNew(rangeNew,rangeOrig) = adj(iOverlapping,:);
adjNew(rangeOrig,rangeNew) = adj(:,iOverlapping);
adjNew(rangeNew,rangeNew) = adj(iOverlapping,iOverlapping);

% Fix the node labels:
commLabelsNew = zeros(newNodes,1);
for i = 1:numNodes
    % If overlapping, take the first label
    theFirstModule = find(commLabels(i,:) > 0,1,'first');
    if ~isempty(theFirstModule)
        commLabelsNew(i) = theFirstModule;
    end
end

% For overlapping, now add the second label
for i = 1:numOverlapping
    theModules = find(commLabels(iOverlapping(i),:) > 0);
    commLabelsNew(rangeNew(i)) = theModules(2);
end