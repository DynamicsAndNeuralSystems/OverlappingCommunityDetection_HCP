function done = generate_synthetic_networks(number, output_dir)
% Function to generate synthetic networks
%
% INPUT:
% number: Number of networks to generate
% Aditi Jha, 29/12/2020
%--------------------------------------------------------

cd weighted_networks

% Check if user supplied output directory
if nargin < 2
    output_dir = getcwd();
end

mkdir(sprintf('%s/networks', output_dir));
mkdir(sprintf('%s/communities', output_dir));

for i = 1:number
   
    % Define new file paths
    newpath_network = sprintf('%s/networks/network%s.dat', output_dir, int2str(i));
    newpath_comm = sprintf('%s/communities/community%s.dat', output_dir, int2str(i));

    % Check if file exists
    if isfile(newpath_network) && isfile(newpath_comm)
        continue
    end
    
    % number of nodes
    N = 180;

    % average degree
    k = 29;

    % maximum degree
    maxk = 102;

    % power law exponent between degree and strength
    beta = 1;

    % power law exponent for degree distribution
    t1 = 2;

    % power law exponent for community sizes
    t2 = randi([2,4],1);

    % Overlapping nodes
    on = randi([0.1*N, 0.15*N],1);

    % membership of overlapping nodes
    om = randi([2,3],1);

    % Mixing parameters
    mut = randi([2,4],1)*0.1;
    muw = randi([2,4],1)*0.1;

    gennet = ['./benchmark -N ', int2str(N),' -k ',int2str(k), ' -maxk ', int2str(maxk), ...
                 ' -beta ', int2str(beta), ' -t1 ', int2str(t1), ' -t2 ', int2str(t2), ...
                 ' -on ', int2str(on), ' -om ', int2str(om), ' -muw ', num2str(muw), ' -mut ', num2str(mut)];
    system(gennet);

    % Move generated files to designated output directory
    system(sprintf('mv network.dat %s', newpath_network));
    system(sprintf('mv community.dat %s', newpath_comm));
end

done = true;
cd('../')

end
