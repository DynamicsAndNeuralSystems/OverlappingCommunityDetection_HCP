function [adjNew, commLabelsNew, iOverlappingNew, nodeListTracking] = constructOverlappingDuplicates(adj, iOverlapping, commLabels)

numNodes = length(adj);
numOverlapping = length(iOverlapping);

% Determine how many new nodes we need based on the number of communities per overlapping node
numDuplicates = sum(sum(commLabels(iOverlapping,:) > 0, 2) - 1);
newNodes = numNodes + numDuplicates;

% Initialize new adjacency matrix
adjNew = zeros(newNodes, newNodes);
adjNew(1:numNodes, 1:numNodes) = adj;

% Track node indices from the original matrix
nodeListTracking = zeros(newNodes, 1);
nodeListTracking(1:numNodes) = 1:numNodes;

% Initialize new community labels
commLabelsNew = zeros(newNodes, 1);
commLabelsNew(1:numNodes) = max(commLabels, [], 2); % Assign first detected label

% Create new indices for duplicate nodes
rangeNewStart = numNodes + 1;
newIndex = rangeNewStart;

iOverlappingNew = false(newNodes, 1);
iOverlappingNew(iOverlapping) = true;

% Map original overlapping nodes to their new duplicates
duplicateMapping = containers.Map('KeyType', 'int32', 'ValueType', 'any');

for i = 1:numOverlapping
    nodeIdx = iOverlapping(i);
    theModules = find(commLabels(nodeIdx, :) > 0);
    numExtraCopies = length(theModules) - 1;

    duplicateIndices = [];
    
    for j = 1:numExtraCopies
        % Assign new index
        nodeListTracking(newIndex) = nodeIdx;
        iOverlappingNew(newIndex) = true;
        
        % Duplicate adjacency information
        adjNew(newIndex, 1:numNodes) = adj(nodeIdx, :);
        adjNew(1:numNodes, newIndex) = adj(:, nodeIdx);
        
        % Assign community label
        commLabelsNew(newIndex) = theModules(j + 1);
        
        % Store duplicate indices for connectivity updates
        duplicateIndices = [duplicateIndices, newIndex];
        
        % Increment index
        newIndex = newIndex + 1;
    end
    
    % Store the mapping for the duplicated nodes
    duplicateMapping(nodeIdx) = duplicateIndices;
end

% Ensure new duplicate nodes are correctly linked
for i = 1:numOverlapping
    nodeIdx = iOverlapping(i);
    if isKey(duplicateMapping, nodeIdx)
        duplicateIndices = duplicateMapping(nodeIdx);
        
        % Connect duplicates to each other
        for j = 1:length(duplicateIndices)
            for k = j+1:length(duplicateIndices)
                adjNew(duplicateIndices(j), duplicateIndices(k)) = adj(nodeIdx, nodeIdx);
                adjNew(duplicateIndices(k), duplicateIndices(j)) = adj(nodeIdx, nodeIdx);
            end
        end
    end
end

end
