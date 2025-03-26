%-------------------------------------------------------------------------------
% Add helper_funcs to path
addpath('helper_funcs');

%-------------------------------------------------------------------------------
% Load processed right-hemisphere HCP structural connectivity data
% Roberts method (0.15) on iFOD
load('../data/HCP_Connectome/RH.mat','RH');

%-------------------------------------------------------------------------------
f = figure('color','w');

% Connectome visualisation
subplot(2,3,[1,2,4,5]); axis('square')
PlotConnectome(RH,true,false);

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
