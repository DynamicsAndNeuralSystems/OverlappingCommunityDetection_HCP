function [adjNew,commLabelsNew,iOverlappingNew] = constructOverlappingDuplicates(adj,iOverlapping,commLabels)

numNodes = length(adj);
numOverlapping = length(iOverlapping);
newNodes = numNodes + numOverlapping;

% Duplicate the main adj
adjNew = zeros(newNodes,newNodes);
adjNew(1:numNodes,1:numNodes) = adj;

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

%-------------------------------------------------------------------------------
% Proxy label for overlapping in new indexing
iOverlappingNew = false(newNodes,1);
iOverlappingNew(iOverlapping) = true;
iOverlappingNew(rangeNew) = true;

% RHnew(1:180,1:180) = RH;
% RHnew(181:191,1:180) = RH(overlapping,:);
% RHnew(1:180,181:191) = RH(:,overlapping);
% RHnew(181:191,181:191) = RH(overlapping,overlapping);

end
