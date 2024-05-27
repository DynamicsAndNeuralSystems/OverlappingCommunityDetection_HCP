% A comparison in runtime for Louvain and OSLOM

rng('default')
storedBenchmark = fullfile(GiveMeFile('OCDA_toolbox'),'Computation','Conversions','networks');
networkDataFile = fullfile(storedBenchmark,'network56.dat');
nodeLabelDataFile = fullfile(storedBenchmark,'community56.dat');


load(networkDataFile,'network56')
% cd(GiveMeFile('OCDA_toolbox'));

DirList = network56;
numNodes = max(max(DirList(:,1:2))); % Calculates the number of nodes in the system

Mat = Direct2Matrix(DirList,numNodes); % Calls function to make matrix input
Undir = Mat2Undir(Mat); % Calls function to make undirected input


REPEAT = 10;
times = zeros(REPEAT,2);

for i = 1:REPEAT
    tic;
    evalc('community_louvain(Mat)');
    times(i,1) = toc;
end

cd(GiveMeFile('OCDA_toolbox'));
for i = 1:REPEAT
    tic;
    evalc('call_OSLOM(Undir, numNodes, 100, 0.1)');
    times(i,2) = toc;
end

sprintf("Louvain Average Time over %d runs = %f seconds with std dev %f",REPEAT,mean(times(:,1)),std(times(:,1)))
sprintf("OSLOM Average Time over %d runs = %f seconds with std dev %f",REPEAT,mean(times(:,2)),std(times(:,2)))