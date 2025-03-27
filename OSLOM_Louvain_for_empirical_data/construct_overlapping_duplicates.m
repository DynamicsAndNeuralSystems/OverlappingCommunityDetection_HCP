function [connectivity_mat_with_repeats, node_indices_with_repeats] = annie_overlapping_duplicates(connectivity_mat, indices_to_repeat, num_of_repeats, node_indices_original)
        % Initialize variables
    num_original_nodes = size(connectivity_mat, 1);
    total_new_nodes = num_original_nodes + sum(num_of_repeats);
    
    % Create an identity mapping of nodes
    node_list = 1:num_original_nodes;
    
    % Prepare new connectivity matrix
    connectivity_mat_with_repeats = zeros(total_new_nodes, total_new_nodes);
    
    % Initialize an empty index map for new positions
    new_order = [];
    
    % Iterate through original nodes and insert duplicates immediately after
    new_idx = 1;
    for original_idx = 1:num_original_nodes
        % Always keep the original node
        new_order = [new_order, original_idx];
        
        % Check if this node should be repeated
        repeat_idx = find(indices_to_repeat == original_idx, 1);
        if ~isempty(repeat_idx)
            num_repeats = num_of_repeats(repeat_idx);
            
            % Insert duplicates immediately after
            for j = 1:num_repeats
                new_order = [new_order, original_idx];
            end
        end
    end
    
    % Reconstruct the connectivity matrix based on new order
    connectivity_mat_with_repeats = connectivity_mat(new_order, new_order);

    % Initialize new vector for node indices
    num_original = length(node_indices_original);
    total_new_length = num_original + sum(num_of_repeats);
    node_indices_with_repeats = zeros(1, total_new_length);
    
    % Build new order with repeated values
    new_order = [];
    
    for i = 1:num_original
        % Always keep the original value
        new_order = [new_order, i];
    
        % Check if this index should be repeated
        repeat_idx = find(indices_to_repeat == i, 1);
        if ~isempty(repeat_idx)
            num_repeats = num_of_repeats(repeat_idx);
            
            % Insert duplicates immediately after
            new_order = [new_order, repmat(i, 1, num_repeats)];
        end
    end
    
    % Construct the new vector based on the expanded order
    node_indices_with_repeats = node_indices_original(new_order);
end
