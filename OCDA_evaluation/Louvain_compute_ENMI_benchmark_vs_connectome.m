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
num_networks = 1000;

% Use a seed of 127 for Louvain
Louvain_seed=127;

% Use default of gamma=1
gamma=1.0;

%%
ENMI_res = load(sprintf('%s/ENMI_Results/all_benchmark_OCDA_ENMI.mat', network_output_dir)).all_ENMI_res;
overlap_sensitivity_res = load(sprintf('%s/ENMI_Results/all_benchmark_OCDA_sensitivity.mat', network_output_dir)).all_overlap_sensitivity_res;
overlap_specificity_res = load(sprintf('%s/ENMI_Results/all_benchmark_OCDA_specificity.mat', network_output_dir)).all_overlap_specificity_res;

%%
for i = 1:num_networks
    % try
        % 
        % Load each synthetic networks previously generated
        benchpath = sprintf('%s/networks/network%s.dat', network_output_dir, int2str(i));
        benchmark = load(benchpath);
        commpath = sprintf('%s/communities/community%s.dat', network_output_dir, int2str(i));
    
        % Load the ground-truth communities and binarize for overlap
        ground_truth_communities = process_Benchmark(commpath, numNodes);
        ground_truth_overlapping = sum(ground_truth_communities ~= 0, 2) >= 2;
    
        % Run OCDAs: networkAdj,methodList,isBenchmark,benchmarkFileName,sourceCodePath, outputPath
        methods_list = {};
        [Final, computation_times] = OLCD_Compute(benchmark, methods_list, 1, ...
            commpath, sourceCodePath, sprintf('%s/ENMI_Results', network_output_dir));
    
        % Extracting all method names
        fields = fieldnames(Final);
        fields = fields(4:end);
        Methods = fields;
    
        % Also compute for Louvain clustering
        rng(Louvain_seed,'twister');
        this_network = Final.Network;
    
        % Widen resulting community matrix
        Louvain_community_affils = community_louvain(this_network, gamma);
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
        Louvain_overlapping = sum(Louvain_matrix_res ~= 0, 2) >= 2;
    
        % Compute ENMI for Louvain
        this_Louvain_ENMI = ENMI_calc(Final.Benchmark.Result, Louvain_matrix_res);
        ENMI_res(23,i) = this_Louvain_ENMI;

        % Compute the sensitivity/specificity in detecting overlaps
        TP = sum((ground_truth_overlapping == 1) & (Louvain_overlapping == 1));
        TN = sum((ground_truth_overlapping == 0) & (Louvain_overlapping == 0));
        FP = sum((ground_truth_overlapping == 0) & (Louvain_overlapping == 1));
        FN = sum((ground_truth_overlapping == 1) & (Louvain_overlapping == 0));
        
        % Sensitivity: True Positive Rate / Recall
        Louvain_sensitivity = TP / (TP + FN);
        overlap_sensitivity_res(23,i) = Louvain_sensitivity;

        % Specificity: True Negative Rate / Precision
        Louvain_specificity = TN / (TN + FP);
        overlap_specificity_res(23,i) = Louvain_specificity;
    % end
end

%% Save results

all_ENMI_res = ENMI_res;
all_overlap_specificity_res = overlap_specificity_res;
all_overlap_sensitivity_res = overlap_sensitivity_res;

% Save final results
save(sprintf('%s/ENMI_Results/all_benchmark_OCDA_ENMI.mat', network_output_dir),'all_ENMI_res');
save(sprintf('%s/ENMI_Results/all_benchmark_OCDA_specificity.mat', network_output_dir),'all_overlap_specificity_res');
save(sprintf('%s/ENMI_Results/all_benchmark_OCDA_sensitivity.mat', network_output_dir),'all_overlap_sensitivity_res');