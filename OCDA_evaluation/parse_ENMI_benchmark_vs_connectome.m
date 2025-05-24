% Add computation code to path
OCDA_evaluation_path = sprintf('%s/OCDA_evaluation', fileparts(pwd));
addpath(genpath(fullfile(OCDA_evaluation_path, 'Computation')));
addpath(genpath(fullfile(OCDA_evaluation_path, 'SourceCode')));
addpath(genpath(fullfile(fileparts(pwd), 'Peripheral')));
sourceCodePath = fullfile(OCDA_evaluation_path, 'SourceCode');

network_output_dir = [fileparts(pwd), '/data/'];
num_networks = 1000;
num_methods = 23;

all_ENMI_res = zeros(num_methods, num_networks);
all_overlap_specificity_res = zeros(num_methods, num_networks);
all_overlap_sensitivity_res = zeros(num_methods, num_networks);
all_OCDA_computation_times = zeros(num_methods, num_networks);

for i = 1:num_networks
    if isfile(sprintf('%s/ENMI_Results/individual_networks/computation_time_network%g.mat', network_output_dir, i))
        % Load network ENMI results
        ENMI_res = load(sprintf('%s/ENMI_Results/individual_networks/ENMI_network%g.mat', network_output_dir, i)).ENMI_res;
       
        overlap_specificity_res = load(sprintf('%s/ENMI_Results/individual_networks/specificity_network%g.mat', network_output_dir, i)).overlap_specificity_res;
        overlap_sensitivity_res = load(sprintf('%s/ENMI_Results/individual_networks/sensitivity_network%g.mat', network_output_dir, i)).overlap_sensitivity_res;

        OCDA_computation_times = load(sprintf('%s/ENMI_Results/individual_networks/computation_time_network%g.mat', network_output_dir, i)).this_comp_time_array;
        OCDA_computation_times = reshape(OCDA_computation_times, 23, 1);

        % Add results to corresponding matrices
        all_ENMI_res(:, i) = ENMI_res;
        all_overlap_specificity_res(:, i) = overlap_specificity_res;
        all_overlap_sensitivity_res(:, i) = overlap_sensitivity_res;
        all_OCDA_computation_times(:, i) = OCDA_computation_times;

    end
end

% Save final results
save(sprintf('%s/ENMI_Results/all_benchmark_OCDA_ENMI.mat', network_output_dir),'all_ENMI_res');
save(sprintf('%s/ENMI_Results/all_benchmark_OCDA_specificity.mat', network_output_dir),'all_overlap_specificity_res');
save(sprintf('%s/ENMI_Results/all_benchmark_OCDA_sensitivity.mat', network_output_dir),'all_overlap_sensitivity_res');
save(sprintf('%s/ENMI_Results/all_OCDA_computation_times.mat', network_output_dir),'all_OCDA_computation_times');