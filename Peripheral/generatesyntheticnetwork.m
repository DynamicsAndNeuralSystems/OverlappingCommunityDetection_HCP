function done = generatesyntheticnetwork(number)
% Function to generate synthetic networks
%
% INPUT:
% number: Number of networks to generate
% Aditi Jha, 29/12/2020
%--------------------------------------------------------

cd weighted_networks
mkdir networks
mkdir communities

for i = 1:number
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
    newpath_net = ['mv network.dat networks/network', int2str(i), '.dat'];
    newpath_comm = ['mv community.dat communities/community', int2str(i), '.dat'];
    system(newpath_net)
    system(newpath_comm)
end

done = true;
cd('../')

end
