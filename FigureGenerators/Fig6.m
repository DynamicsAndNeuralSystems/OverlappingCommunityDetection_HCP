% load RH connectome
load('RH.mat');

% Load Louvain community labels for RH connectome
load('louvaincomm.mat');

%Calculate Participation Coefficient
P = participation_coef(RH,louvaincomm);

% Calculate within-module degree z-scored
z = module_degree_zscore(RH,louvaincomm);

% subplot(1,2,1); axis square;

color1 = [0.66,0.66,0.66];
scatter(P,z,80,color1,'filled');
hold ('on')

% Load OSLOM community affiliation vectors on RH connectome
load('oslomcomm.mat');
iOverlapping = find(sum(oslomcomm > 0, 2) == 2); % list of overlapping nodes obtained from OSLOM

scatter(P(iOverlapping), z(iOverlapping), 80,'filled','r');
xlabel('Participation Coefficient, P');
ylabel('Within-module strength, z');
xlim([0 1]);
legend('Non-Overlapping nodes','Overlapping nodes','Location','southeast');
title('Louvain','Fontsize',15);

%-------------------------------------------------------------------------------
% For OSLOM, replacing each overlapping node by two non-overlapping nodes
% (since every overlapping node belong to two communities)
%-------------------------------------------------------------------------------
[RHnew,oslomRHnew,iOverlappingNew] = constructOverlappingDuplicates(RH,iOverlapping,oslomcomm);

%-------------------------------------------------------------------------------
% Load OSLOM community labels for RHnew connectome
% load('oslomRHnew.mat','oslomRHnew');

% Calculate Participation Coefficient
P = participation_coef(RHnew,oslomRHnew);

% Calculate within-module degree z-scored
z = module_degree_zscore(RHnew,oslomRHnew);

% subplot(1,2,2); axis square;
figure('color','w');
color1 = [0.66,0.66,0.66];
scatter(P,z,80,color1,'filled');
hold('on')

scatter(P(iOverlappingNew), z(iOverlappingNew), 80,'filled','r');
xlabel('Participation Coefficient, P');
ylabel('Within-module strength, z');
xlim([0 1]);
legend('Non-Overlapping nodes','Overlapping nodes','Location','southeast');
title('OSLOM','Fontsize',15);
