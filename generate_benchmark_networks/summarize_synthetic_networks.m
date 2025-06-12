% Generate 1000 synthetic networks using generate_synthetic_networks.m,
% Then compute extended normalized mutual information (ENMI) for
% each of the 22 atatlmethods on each of the 1000 synthetic networks
% Two arguments to generate_synthetic_networks: number of networks to generate and the directory to save the networks
num_networks = 1000;
numNodes = 180;
network_output_dir = [fileparts(pwd), '/data/'];

benchmark_ensemble_summary_stats = zeros(num_networks, 4);

for i = 1:num_networks

    % Load each synthetic networks previously generated
    benchpath = sprintf('%s/networks/network%s.dat', network_output_dir, int2str(i));
    benchmark = load(benchpath);
    commpath = sprintf('%s/communities/community%s.dat', network_output_dir, int2str(i));

    % Load the ground-truth communities and binarize for overlap
    ground_truth_communities = process_Benchmark(commpath, numNodes);

    % Define (1) number of communities, (2) number of overlapping nodes,
    % (3) size of largest community, (4) size of smallest community

    % (1) number of communities
    number_of_communities = size(ground_truth_communities, 2);

    % (2) number of overlapping nodes
    num_overlapping_nodes = sum(sum(ground_truth_communities ~= 0, 2) >= 2);

    % size of (3) largest and (4) smallest community
    community_sizes = sum(ground_truth_communities~=0, 1);
    largest_community_size = max(community_sizes);
    smallest_community_size = min(community_sizes);

    % Add data points to matrix
    benchmark_ensemble_summary_stats(i,:) = [number_of_communities, num_overlapping_nodes, ... 
                                                largest_community_size, smallest_community_size];

end

% Save ensemble summary statistics to a CSV
writematrix(benchmark_ensemble_summary_stats, ...
    sprintf('%s/Benchmark_ensemble_summary_stats.csv', network_output_dir));