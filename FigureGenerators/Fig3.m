%-------------------------
% ------------------------------------------------------
% Load benchmark network data
storedBenchmark = fullfile(GiveMeFile('OCDA_toolbox'),'Computation','Conversions','networks');
networkDataFile = fullfile(storedBenchmark,'network56.dat');
nodeLabelDataFile = fullfile(storedBenchmark,'community56.dat');

% Run OCDAs on this benchmark
load(networkDataFile,'network56')
cd(GiveMeFile('OCDA_toolbox'));
OLCD_Compute(network56, {'OSLOM','Jerry','Shen','NNMF','Infomap','SLPA'}, true, nodeLabelDataFile)
% OLCD_Compute(network56, {'Shen','NNMF','Infomap'}, true, nodeLabelDataFile)

% Visualize performance
% cd ../
% cd Visualization_Module
OLCD_Visualization({'Benchmark', 'OSLOM', 'Shen', 'NNMF', 'Infomap', 'SLPA'})


% ORIGINAL ADITI:
% %load benchmark network data
% generatesyntheticnetwork(1)
% network = load('weighted_networks/networks/network1.dat');

% cd /Users/aditijha/Desktop/communityDetection/
% %Run OCDAs on benchmark
% Computation(network, {'OSLOM', 'Clique', 'NNMF', 'Infomap','SLPA'}, 1, 'weighted_networks/communities/community1');

% % Visualize performance
% % cd ../
% % cd Visualization_Module
% Visualization({'Benchmark', 'OSLOM', 'Shen', 'NNMF', 'Infomap', 'SLPA'})
