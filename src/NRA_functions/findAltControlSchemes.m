function [to_manipulate_all,manipulate_changes_all,len_manipulate,DPs,objectives,model] = ...
    findAltControlSchemes(model,options,ind_bfuse,consnames)

% This functions finds noAlt alternatives for the specific model
[DPs,objectives,model] = findDP_CONTROL_NEW(model,options.noAlt, ind_bfuse,options.timelimit);
for i=1:size(DPs,2)
    DP = DPs(:,i);
    % Find the bfuse vars
    cont = DP(ind_bfuse);
    % whichc enzymes to change?
    to_manipulate = consnames(find(cont<0.1));
    % store ones to manipulate
    to_manipulate_all{i} = to_manipulate;
    % calculate net change F - R
    net_change = DP(find_cell(strcat('F_',consnames),model.varNames))-...
        DP(find_cell(strcat('R_',consnames),model.varNames));
    % Manipulate changes compute and store
    manipulate_changes = net_change(find_cell(to_manipulate,consnames));    
    manipulate_changes_all{i} = manipulate_changes;
    len_manipulate(i) = length(to_manipulate_all{i});
end

end