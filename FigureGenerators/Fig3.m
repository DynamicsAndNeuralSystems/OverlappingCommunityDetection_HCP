%-------------------------------------------------------------------------------
% Load benchmark network data
storedBenchmark = fullfile(GiveMeFile('OCDA_toolbox'),'Computation','Conversions','networks');
networkDataFile = fullfile(storedBenchmark,'network56.dat');
nodeLabelDataFile = fullfile(storedBenchmark,'community56.dat');

% Run OCDAs on this benchmark
load(networkDataFile,'network56')
cd(GiveMeFile('OCDA_toolbox'));
Computation(network56, {'OSLOM','Jerry','Shen','NNMF','Infomap'}, true, nodeLabelDataFile)

% Visualize performance
% cd ../
% cd Visualization_Module
Visualization({'Benchmark', 'OSLOM', 'Shen', 'NNMF', 'Infomap', 'SLPA'})
