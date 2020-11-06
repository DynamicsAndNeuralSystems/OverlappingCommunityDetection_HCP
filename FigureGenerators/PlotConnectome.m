function Plot = PlotConnectome(adjMatrix,doNormalize,doReorder)
% Plot the connectome contained in the adjacency matrix, adjMatrix

if nargin < 2
    doNormalize = false;
end
if nargin < 3
    doReorder = true;
end
%-------------------------------------------------------------------------------

if doNormalize
    adjMatrix = adjMatrix/max(adjMatrix(:));
end
if doReorder
    ix = BF_ClusterReorder(adjMatrix,'euclidean','average');
    adjMatrix = adjMatrix(ix,ix);
end
adjMatrix(adjMatrix == 0) = NaN; % Turns the 0s to NaNs for easy plotting
adjMatrix = flipud(adjMatrix);

%% Figure options
hold('on')
Plot = imagesc(adjMatrix);
axis('image');
set(gca, 'color', [0 0 0], 'CLim', [0 max(adjMatrix(:))]); % Sets background colour to black
set(Plot, 'alphadata',~isnan(adjMatrix)); % Turns all NaN values into transparent colours
c = colorbar('location', 'eastoutside'); % Shows a colour bar
colormap(GiveMeTurbo());

%% Labeling
title('A','FontSize', 15);
XTickLabel = size(adjMatrix, 1)+0.5; % Kickstarts the for loop
set(gcf, 'InvertHardCopy', 'off'); % Fixes white background issue
c.Label.String = "Edge Weight";
c.Label.FontSize = 10;

xlabel('Region')
ylabel('Region')

end
