% -------------------------------------------------------------------------------
%% REPRODUCING FIG 8
%-------------------------------------------------------------------------------

%% LOAD DATA
% Loading Right hemisphere connectome
load('RH.mat')
% Loading community structure as determined by Louvain
load('louvaincomm_Rubinov.mat');
% Loading community structure as determined by OSLOM
load('oslomcomm.mat');
% load('oslomRHnew.mat');

%-------------------------------------------------------------------------------
%% FILTER
%-------------------------------------------------------------------------------
% Keeping all edges that connect to node 63 (which is the target node to visualize based on Fig. 7)
keepNode = RH(63,:) > 0;
keepNode(63) = 1;
RHsub = RH(keepNode, keepNode);
numNodes = size(RHsub,1);

% Filter community assignments by Louvain and OSLOM
lcomm = louvaincomm_Rubinov(keepNode);
ocomm = oslomcomm(keepNode,:);

%-------------------------------------------------------------------------------
%% LIST
%-------------------------------------------------------------------------------
% Creating a list with community IDs since OSLOM output is in the form of
% one-hot matrix
max_comm = size(ocomm,2);
ocomm_list = zeros(numNodes,1);
for i = 1:numNodes
    % Community assignments for ith node
    comms_i = find(ocomm(i,:));
    ncomms_i = size(comms_i,2);
    if ncomms_i > 1
        % Overlapping node between networks 1 + 2:
        if ocomm(i,1) > 0 && ocomm(i,2) > 0
            ocomm_list(i) = max_comm + 1;
        elseif ocomm(i,1) > 0 || ocomm(i,2) > 0
                ocomm_list(i) = max_comm + 2;
        end
    else
        ocomm_list(i) = comms_i(1);
    end
end

nodelist = (1:numNodes);

%-------------------------------------------------------------------------------
%% PLOT
%-------------------------------------------------------------------------------
figure()
% OSLOM
% subplot(1,2,2)
G = graph(RHsub);
p = plot(G);


% Communities 1 and 2 are merged by Louvain
% nodes_overlapping_1_and_2 represent overlapping nodes (between communities 1 and 2)
% nodes_overlapping_other represent overlapping nodes that involve communities 1 or 2 and a different community

nodes1 = nodelist(ocomm_list==1);
nodes2 = nodelist(ocomm_list==2);
nodes_overlapping_1_and_2 = nodelist(ocomm_list==7);
nodes_overlapping_other = nodelist(ocomm_list==8);

highlight(p,nodelist, 'NodeColor', [0.9,0.9,0.9]);
highlight(p,nodes1, 'NodeColor', [100/256, 199/256, 232/256]);
highlight(p,nodes2, 'NodeColor', [200/256, 127/256, 14/256]);
highlight(p,nodes_overlapping_1_and_2, 'NodeColor', [128/255, 0, 0]);
highlight(p,nodes_overlapping_other, 'NodeColor', [0, 0, 0]);

p.MarkerSize = 20;
p.EdgeColor = [0.7, 0.7, 0.7];
p.NodeFontSize = 0.5;

% Swap the x and y axis data
xData = get(p, 'XData');
yData = get(p, 'YData');

% Assign the swapped data
set(p, 'XData', yData);
set(p, 'YData', xData);

colormap(jet);
title('OSLOM');

% LOUVAIN
% subplot(1,2,1)
figure()
G = graph(RHsub);
p = plot(G);

p.MarkerSize = 20;
p.EdgeColor = [0.7, 0.7, 0.7];
% Highlighting the nodes belong to community 2 which is split into several
% communities by OSLOM
nodes2 = nodelist(lcomm==2);
highlight(p,nodelist, 'NodeColor', [0.9,0.9,0.9]);
highlight(p,nodes2, 'NodeColor', [100/256, 199/256, 232/256]);
% highlight(p,nodes7, 'NodeColor', [128/255, 0, 0]);
p.NodeFontSize = 0.5;
colormap(jet);
title('Louvain');

% Swap the x and y axis data
xData = get(p, 'XData');
yData = get(p, 'YData');

% Assign the swapped data
set(p, 'XData', yData);
set(p, 'YData', xData);



numnodes = size(keepNode);
lookupmatrix = -1 * ones(2, numnodes(2));
lookupmatrix(1,:) = keepNode;

currentidx = 1;
for i = 1:numnodes(2)
    if keepNode(i)==1
        lookupmatrix(2,i) = currentidx;
        currentidx = currentidx +1 ;
    end
end


writematrix(lookupmatrix',"../Results/lookupmatrix.csv")