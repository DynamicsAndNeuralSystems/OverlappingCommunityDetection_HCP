% Generate 1000 synthetic networks using generate_synthetic_networks.m,
% Then compute extended normalized mutual information (ENMI) for
% each of the 22 atatlmethods on each of the 1000 synthetic networks
% Two arguments to generate_synthetic_networks: number of networks to generate and the directory to save the networks
num_networks = 1000;
network_output_dir = [fileparts(pwd), '/data/'];

% Store the starting time 
start_time = tic;
generate_synthetic_networks(num_networks, network_output_dir);
% Display the time taken to generate the networks
fprintf('Time taken to generate %d networks: %.2f seconds\n', num_networks, toc(start_time));
duration = toc(start_time);
writematrix(duration, sprintf('%s/generate_networks_time_to_compute.txt', network_output_dir));
