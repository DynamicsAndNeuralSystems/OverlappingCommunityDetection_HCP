function [adjNew,commLabelsNew,iOverlappingNew] = create_duplicated_adj(RH_in, RH_out, oslom_labels, github_path)
    % Add github path to path
    addpath(genpath(github_path));

    % Load in RH_in
    load(RH_in);

    % Load in OSLOM labels
    load(oslom_mat);

    iOverlapping = find(sum(oslomcomm > 0, 2)==2); % list of overlapping nodes obtained from OSLOM

    [RH_new,oslomRHnew,iOverlappingNew,nodeListTracking, commLabelsNew] = constructOverlappingDuplicates(RH,iOverlapping,oslomcomm);

    % Save RH_new to .mat
    save(RH_out, "oslomRHnew", "iOverlappingNew", "RH_new","nodeListTracking", "commLabelsNew");
end

