%-------------------------------------------------------------------------------
% Load processed right-hemisphere HCP structural connectivity data
% Roberts method (0.15) on iFOD
load(fullfile('Data','RH.mat'),'RH')

%-------------------------------------------------------------------------------
f = figure('color','w');

% Connectome visualisation
subplot(2,3,[1,2,4,5]); axis('square')
PlotConnectome(RH,false,true);

%-------------------------------------------------------------------------------
% Degree vs strength plot
subplot(2,3,3);
binRH = (RH > 0);
degree = sum(binRH,1);
strength = sum(RH,1);
plot(degree,strength,'ok','MarkerFaceColor','k','MarkerSize',3)
xlim([0,max(degree)])
title('B','FontSize', 15);
xlabel('Degree, k')
ylabel('Strength, s')

%-------------------------------------------------------------------------------
% Degree Histogram
subplot(2,3,6);
axis('square')
histogram(degree,20,'FaceColor',0.5*ones(1,3),'EdgeColor','k')
xlim([0,max(degree)])
title('C','FontSize', 15);
xlabel('Degree, k')
ylabel('Frequency')

f.Position = [1000        1002         648         336];

%-------------------------------------------------------------------------------
% Powerlaw?:
% if any(degree==0)
%     numLone = sum(degree==0);
%     warning('%u nodes have no connections--excluding',numLone)
%     degree = degree(degree>0);
% end

% [binCenters,Nnorm] = binLogLog(degree,10);
% f = figure('color','w');
% loglog(binCenters,Nnorm,'ok')
