%load benchmark network data
generatesyntheticnetwork(1)
network = load('weighted_networks/networks/network1.dat');

cd /Users/aditijha/Desktop/communityDetection/
%Run OCDAs on benchmark
Computation(network, {'OSLOM', 'Clique', 'NNMF', 'Infomap','SLPA'}, 1, 'weighted_networks/communities/community1');

%Visualize performance
cd Visualisation_Module
Visualisation({'Benchmark', 'OSLOM', 'Clique', 'NNMF', 'Infomap','SLPA'})

