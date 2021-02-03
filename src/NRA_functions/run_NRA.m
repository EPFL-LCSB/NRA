function runEcoliD1_NRA_PCA_QP(save_dir,model_name,options,NRA_dir)
% Function for running the NRA on different cases.

caseStudy = 1;

switch caseStudy
    case 1
        % Load and concatinate the CCs, kinModel etc
        load([save_dir,'metSampData',model_name,'.mat'],'kinModel')
        load([save_dir,'CCforNRA_',options.objrxnPCA{:},'.mat']);
        controlflux = CCforNRA.aggFCC(:,:,1);
        controlconc = CCforNRA.aggCCC(:,:,1);
        consnames = CCforNRA.consnames;
        fluxnames = CCforNRA.fluxnames;
        % Create the NRA model
        [NRAmodel, ind_bfuse] = prepMCAcoefficientsForNRA(controlflux, controlconc, consnames,kinModel, options, fluxnames);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fix glucose splitting ratios based on reference data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(options,'fixGLCptsppRatio') && options.fixGLCptsppRatio == 1
    [num_constr, ~] = size(NRAmodel.A);
    NRAmodel.constraintNames{num_constr+1,1} = 'GLCratio';
    NRAmodel.constraintType{num_constr+1,1} = '=';
    NRAmodel.rhs(num_constr+1) = 0;
    NRAmodel.A(num_constr+1,find_cell('RXN_GLCptspp',NRAmodel.varNames)) = 1;
    NRAmodel.A(num_constr+1,find_cell('RXN_GLCtex',NRAmodel.varNames)) = -1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust bounds for fluxes/concentrations for changes that we know or want
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NRAmodel.var_lb(find_cell(options.objvars, NRAmodel.varNames)) = options.objvlb;
NRAmodel.var_ub(find_cell(options.objvars, NRAmodel.varNames)) = options.objvub;
sol = solveTFAmodelCplex(NRAmodel,300);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the objective function and problem type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(options.setobjective,'Rxn')
    NRAmodel.f = zeros(length(NRAmodel.varNames),1);
    NRAmodel.f(find_cell(options.objrxn,NRAmodel.varNames)) = 1;
    NRAmodel.objtype = options.objtype;
    % Check to see if there is a solution with a small time limit
    sol = solveTFAmodelCplex(NRAmodel,300);
end

if ~isempty(sol.val) && strcmp(options.setobjective,'Enzymes')
    ind = getAllVar(NRAmodel,{'BFUSE'});
    manipulate_vars = NRAmodel.varNames(ind(find(sol.x(ind)<0.1)));
    manipulate_vars_sep = [strrep(manipulate_vars,'BFUSE','F') ; strrep(manipulate_vars,'BFUSE','R')];
    sol_manipulate = sol.x(find_cell(manipulate_vars_sep,NRAmodel.varNames));
    sol_manipulate_idx = find(abs(sol_manipulate)>1.e-9);
    to_manipulate = [manipulate_vars_sep(sol_manipulate_idx) num2cell(sol_manipulate(sol_manipulate_idx))];
    [~,idx] = sort(abs(sol_manipulate(sol_manipulate_idx)),'descend');
    to_manipulate_sorted = to_manipulate(idx,:);
    save([NRA_dir,'to_manipulate_sorted_ref.mat'],'to_manipulate_sorted')
elseif isempty(sol.val)
    % No solution checkpoint
    disp('Infeasible model! Check your bounds!')
    keyboard
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([NRA_dir,'NRAmodel.mat'],'NRAmodel')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enumerate alternative strategies. Note here that it searches for up to
% options.noAlt. Once a no solution alternative is found it will stop. If
% the number of allowed enzymatic changes options.num_enz_perturb is not
% fixed, this can run for very long.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.enumAlternatives == 1 && options.runAllSaturations == 0
    [to_manipulate_all,manipulate_changes_all,len_manipulate,DPs,objectives] = ...
        findAltControlSchemes(NRAmodel,options,ind_bfuse,consnames);
    DataNRA.to_manipulate_all = to_manipulate_all;
    DataNRA.manipulate_changes_all = manipulate_changes_all;
    DataNRA.len_manipulate = len_manipulate;
    DataNRA.DPs = DPs;
    DataNRA.objectives = objectives;
    save([NRA_dir,'DataNRA.mat'],'DataNRA')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Go through all the (19 default) saturations selected from PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.runAllSaturations ==1
    modelCount = 0;
    sol_count = 0;
    has_sol = zeros(size(CCforNRA.aggFCC,3),1);
    all_sols = cell(size(CCforNRA.aggFCC,3),1);
    if options.enumAlternatives == 1
        agg_sols = cell(1,size(CCforNRA.aggFCC,3));
    end
    % Check each selected CC
    for k = 1:size(CCforNRA.aggFCC,3) % num CCs
        
        controlflux = CCforNRA.aggFCC(:,:,k);
        controlconc = CCforNRA.aggCCC(:,:,k);
        
        % Update in the A matrix (aka the N in VH terms) the control
        % coeffcients by the new ones.
        tempMat = [controlconc -controlconc; controlflux -controlflux];
        numConAndFlux = size(controlconc,1)+size(controlflux,1);
        
        % SOLVE THE MILP PROBLEM
        load([NRA_dir,'NRAmodel.mat'])
        NRAmodel.A(1:numConAndFlux, numConAndFlux+1:numConAndFlux+2*numel(consnames)) = tempMat;
        NRAmodel.var_lb(find_cell(options.objvars, NRAmodel.varNames)) = options.objvlb;
        NRAmodel.var_ub(find_cell(options.objvars, NRAmodel.varNames)) = options.objvub;
        
        % Check the objective function and problem type.
        if strcmp(options.setobjective,'Rxn')
            NRAmodel.f = zeros(length(NRAmodel.varNames),1);
            NRAmodel.f(find_cell(options.objrxn,NRAmodel.varNames)) = 1;
            NRAmodel.objtype = options.objtype;
            sol = solveTFAmodelCplex(NRAmodel,options.timelimit);
        end
        
        if ~isempty(sol.val)
            sol_count = sol_count+1;
            has_sol(1,k) = 1;
            all_sols{1,k} = sol;
            ind = getAllVar(NRAmodel,{'BFUSE'});
            manipulate_vars = NRAmodel.varNames(ind(find(sol.x(ind)==0)));
            manipulate_vars_sep = [strrep(manipulate_vars,'BFUSE','F') ; strrep(manipulate_vars,'BFUSE','R')];
            sol_manipulate = sol.x(find_cell(manipulate_vars_sep,NRAmodel.varNames));
            sol_manipulate_idx = find(abs(sol_manipulate)>0);
            to_manipulate = [manipulate_vars_sep(sol_manipulate_idx) num2cell(sol_manipulate(sol_manipulate_idx))];
            [~,idx] = sort(abs(sol_manipulate(sol_manipulate_idx)),'descend');
            to_manipulate_sorted_all.(['case',num2str(k)]) = to_manipulate(idx,:);
            if options.enumAlternatives == 1
                % Enumerate alternative strategies. Note here that it searches for up to
                % options.noAlt. Once a no solution alternative is found it will stop. If
                % the number of allowed enzymatic changes options.num_enz_perturb is not
                % fixed, this can run for very long.
                [to_manipulate_all,manipulate_changes_all,len_manipulate,DPs] = ...
                    findAltControlSchemes(NRAmodel,options,ind_bfuse,consnames);
                DataNRA.to_manipulate_all = to_manipulate_all;
                DataNRA.manipulate_changes_all = manipulate_changes_all;
                DataNRA.len_manipulate = len_manipulate;
                DataNRA.DPs = DPs;
                agg_sols{1,k} = DataNRA;
            end
        end
    end
    % Store all the solutions from NRA runs for this batch
    if options.enumAlternatives == 0
        save([NRA_dir,'DataNRAbatch.mat'],'has_sol','all_sols','to_manipulate_sorted_all')
    else
        save([NRA_dir,'DataNRAbatch.mat'],'has_sol','all_sols','to_manipulate_sorted_all','agg_sols')
    end
end

end
