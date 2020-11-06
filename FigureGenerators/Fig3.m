%-------------------------------------------------------------------------------
% Load benchmark network data
storedBenchmark = fullfile(GiveMeFile('OCDA_toolbox'),'Computation','Conversions','networks');
networkDataFile = fullfile(storedBenchmark,'network56.dat');
nodeLabelDataFile = fullfile(storedBenchmark,'community56.dat');

% Run OCDAs on this benchmark
load(networkDataFile,'network56')
cd(GiveMeFile('OCDA_toolbox'));
Computation(network56, {'OSLOM','Jerry','Shen','NNMF','Infomap'}, true, nodeLabelDataFile)

% Visualise performance
% cd ../
% cd Visualisation_Module
Visualisation({'Benchmark', 'OSLOM', 'Shen', 'NNMF', 'Infomap','SLPA'})
