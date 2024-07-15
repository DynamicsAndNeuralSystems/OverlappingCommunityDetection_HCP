% load RH connectome
load('RH.mat');

%load Louvain community labels for RH connectome
load('louvaincomm_Rubinov.mat');
% This was generated using fast-consensus code from Tandon et. al. (the
% code is no more publically available)
% The other option is to use Mika Rubinov's code:
% peripheral/community_louvain(RH, 1);

%Calculate Participation Coefficient
P = participation_coef(RH,louvaincomm_Rubinov);
louvain_P = P;

% Calculate within-module degree z-scored
z = module_degree_zscore(RH, louvaincomm_Rubinov);
louvain_z = z;

% subplot(1,2,1); axis square;

color1 = [0.66,0.66,0.66];
scatter(P,z,80,color1,'filled');
hold('on')

%load OSLOM community affiliation vectors on RH connectome
%This can be generated using:
% Computation(RH, {'OSLOM'}, 0);
load('oslomcomm.mat');
iOverlapping = find(sum(oslomcomm > 0, 2)==2); % list of overlapping nodes obtained from OSLOM

scatter(P(iOverlapping), z(iOverlapping), 80,'filled','r');
xlabel('Participation Coefficient, P');
ylabel('Within-module strength, z');
xlim([0 1]);
legend('Non-Overlapping nodes','Overlapping nodes','Location','southeast');
title('Louvain','Fontsize',15);

% For OSLOM, replacing each overlapping node by two non-overlapping nodes
% (since every overlapping node belong to two communities)
[RHnew,oslomRHnew,iOverlappingNew] = constructOverlappingDuplicates(RH,iOverlapping,oslomcomm);

%Calculate Participation Coefficient
P = participation_coef(RHnew,oslomRHnew);
oslom_P = P;

% Calculate within-module degree z-scored
z = module_degree_zscore(RHnew, oslomRHnew);
oslom_z = z;

% subplot(1,2,2); axis square;
figure('color','w');
color1 = [0.66,0.66,0.66];
scatter(P,z,80,color1,'filled');
hold('on')

z(iOverlappingNew)

scatter(P(iOverlappingNew), z(iOverlappingNew), 80,'filled','r');
xlabel('Participation Coefficient, P');
ylabel('Within-module strength, z');
xlim([0 1]);
legend('Non-Overlapping nodes','Overlapping nodes','Location','southeast');
title('OSLOM', 'Fontsize',15);
