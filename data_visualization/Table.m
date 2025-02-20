% -------------------------------------------------------------------------------
%% TABLE OF NODES AND THEIR MODULES
%-------------------------------------------------------------------------------

%% LOAD DATA
load('RH.mat')
load('louvaincomm_Rubinov.mat');
load('oslomcomm.mat');

node_labels = readtable("brain_area_labels.csv",'Format','%d%s','Delimiter',',');
node_labels = node_labels(:,2); % Cut off the irrelevant index, use the position instead.


modulenames = ["Occipito-temporal", "Somatomotor", "Ventral Prefrontal",
    "Auditory", "Dorsal Prefrontal", "Brainy McBrainFace"];

output = cell(size(oslomcomm,1),3);


for node = 1:size(oslomcomm,1)
    modules = find(oslomcomm(node,:)); % Get the relevant modules by index which the node is associated with
    if isempty(modules)==1
        output{node,2} = "Unassigned";
        output{node,3} = "Unassigned";
    elseif isscalar(modules)
        output{node,2} = modulenames(modules);
        output{node,3} = "Not";
    else
        output{node,2} = strjoin(modulenames(modules)," + ");
        output{node,3} = "Overlapping";
    end

    output{node,1} = node_labels{node,'Var2'}{1};
end
writecell(output,"data/table_of_modules.csv")