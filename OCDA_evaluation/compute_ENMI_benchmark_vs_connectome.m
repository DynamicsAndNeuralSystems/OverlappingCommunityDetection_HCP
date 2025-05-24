% Add computation code to path
OCDA_evaluation_path = sprintf('%s/OCDA_evaluation', fileparts(pwd));
addpath(genpath(fullfile(OCDA_evaluation_path, 'Computation')));
addpath(genpath(fullfile(OCDA_evaluation_path, 'SourceCode')));
addpath(genpath(fullfile(fileparts(pwd), 'Peripheral')));
sourceCodePath = fullfile(OCDA_evaluation_path, 'SourceCode');

network_output_dir = [fileparts(pwd), '/data/'];

% Matrix to store ENMI values, with an extra column for Louvain
numNodes = 180;
numMethods = 23;
num_networks = 1000

% Use a seed of 127 for Louvain
Louvain_seed=127;

% Use default of gamma=1
gamma=1.0;

%%
for i = 1:num_networks
    if ~isfile(sprintf('%s/ENMI_Results/individual_networks/computation_time_network%g.mat', network_output_dir, i))
        try
            ENMI_res = zeros(numMethods);
            overlap_sensitivity_res = zeros(numMethods);
            overlap_specificity_res = zeros(numMethods);
            computation_time_res = zeros(numMethods);
        
            % 
            % Load each synthetic networks previously generated
            benchpath = sprintf('%s/networks/network%s.dat', network_output_dir, int2str(i));
            benchmark = load(benchpath);
            commpath = sprintf('%s/communities/community%s.dat', network_output_dir, int2str(i));
        
            % Load the ground-truth communities and binarize for overlap
            ground_truth_communities = process_Benchmark(commpath, numNodes);
            ground_truth_overlapping = sum(ground_truth_communities ~= 0, 2) >= 2;
        
            % Run OCDAs: networkAdj,methodList,isBenchmark,benchmarkFileName,sourceCodePath, outputPath
            methods_list = {'Infomap', 'SLPA', 'OSLOM', 'Clique','NNMF'};
            [Final, computation_times] = OLCD_Compute(benchmark, methods_list, 1, ...
                commpath, sourceCodePath, sprintf('%s/ENMI_Results', network_output_dir));
        
            % Extracting all method names
            fields = fieldnames(Final);
            fields = fields(4:end);
            Methods = fields;
        
            % Compute ENMI, sensitivity, and specificity for each method
            for j=1:size(Methods,1)
                % Load the community predictions with the given method
                this_OCDA_comm_res = Final.(Methods{j}).Result;
                this_OCDA_overlapping = sum(this_OCDA_comm_res ~= 0, 2) >= 2;
        
                % Compute ENMI
                ENMI_res(j) = ENMI_calc(Final.Benchmark.Result, this_OCDA_comm_res);
                
                % Compute the sensitivity/specificity in detecting overlaps
                TP = sum((ground_truth_overlapping == 1) & (this_OCDA_overlapping == 1));
                TN = sum((ground_truth_overlapping == 0) & (this_OCDA_overlapping == 0));
                FP = sum((ground_truth_overlapping == 0) & (this_OCDA_overlapping == 1));
                FN = sum((ground_truth_overlapping == 1) & (this_OCDA_overlapping == 0));
                
                % Sensitivity: True Positive Rate / Recall
                sensitivity = TP / (TP + FN);
                overlap_sensitivity_res(j) = sensitivity;
        
                % Specificity: True Negative Rate / Precision
                specificity = TN / (TN + FP);
                overlap_specificity_res(j) = specificity;
                    
            end
        
            % Also compute for Louvain clustering
            rng(Louvain_seed,'twister');
            this_network = Final.Network;
            start_time = tic;
            Louvain_community_affils = community_louvain(this_network, gamma);
            Louvain_duration = toc(start_time);
            computation_times.Louvain = Louvain_duration;
        
            % Widen resulting community matrix
            Louvain_community_affils = Louvain_community_affils(:);
            n_nodes = length(Louvain_community_affils);
            communities = unique(Louvain_community_affils);
            n_communities = length(communities);
        
            % Map original labels to consecutive integers
            [~, ~, new_labels] = unique(Louvain_community_affils);
        
            % Create the binary indicator matrix
            Louvain_matrix_res = zeros(n_nodes, n_communities);
            idx = sub2ind(size(Louvain_matrix_res), (1:n_nodes)', new_labels);
            Louvain_matrix_res(idx) = 1;
        
            % Compute ENMI for Louvain
            this_Louvain_ENMI = ENMI_calc(Final.Benchmark.Result, Louvain_matrix_res);
            ENMI_res(23) = this_Louvain_ENMI;

            % Compute the sensitivity/specificity in detecting overlaps
            TP = sum((ground_truth_overlapping == 1) & (Louvain_overlapping == 1));
            TN = sum((ground_truth_overlapping == 0) & (Louvain_overlapping == 0));
            FP = sum((ground_truth_overlapping == 0) & (Louvain_overlapping == 1));
            FN = sum((ground_truth_overlapping == 1) & (Louvain_overlapping == 0));
            
            % Sensitivity: True Positive Rate / Recall
            Louvain_sensitivity = TP / (TP + FN);
            overlap_sensitivity_res(23) = Louvain_sensitivity;
    
            % Specificity: True Negative Rate / Precision
            Louvain_specificity = TN / (TN + FP);
            overlap_specificity_res(23) = Louvain_specificity;
        
            % Save ENMI, sensitivity, and specificity data for this network
            save(sprintf('%s/ENMI_Results/individual_networks/ENMI_network%g.mat', network_output_dir, i),'ENMI_res');
            save(sprintf('%s/ENMI_Results/individual_networks/sensitivity_network%g.mat', network_output_dir, i),'overlap_sensitivity_res');
            save(sprintf('%s/ENMI_Results/individual_networks/specificity_network%g.mat', network_output_dir, i),'overlap_specificity_res');
        
            % Save computation time results 
            this_comp_time_table = struct2table(computation_times);
            this_comp_time_table = removevars(this_comp_time_table, 'Date');
            this_comp_time_array = table2array(this_comp_time_table);
            save(sprintf('%s/ENMI_Results/individual_networks/computation_time_network%g.mat', network_output_dir, i),'this_comp_time_array');
        end
    end
end

