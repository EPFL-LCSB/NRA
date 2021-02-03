function ConcViolations = runConcViolationStudies(save_dir,model_name,options,NRA_dir,ConcViolations)
% Function for runnning NRA Concentration Violation Studies.
% Preset objective values: from log(1.1) to log(2) for increasing
% objectives and from log(0.9) to log(0) for decreasing objectives.

% Load and concatenate the CCs, kinModel etc
load([save_dir,'metSampData',model_name,'.mat'], 'tFBAmodel')
load([save_dir,'metSampData',model_name,'.mat'], 'kinModel')

load([save_dir,'CCforNRA_',options.objrxnPCA{:},'.mat'])

% CCs and their names for preparing NRA model
controlconc = CCforNRA.aggCCC(:,:,1);
controlflux = CCforNRA.aggFCC(:,:,1);
consnames =  CCforNRA.consnames;
fluxnames = CCforNRA.fluxnames;

options.num_enz_perturb = 1; %To make sure the constraint is added
NRAmodel = prepMCAcoefficientsForNRA(controlflux, controlconc, consnames, kinModel, options, fluxnames);

if isfield(options,'fixGLCptsppRatio') && options.fixGLCptsppRatio == 1
    [num_constr, ~] = size(NRAmodel.A);
    NRAmodel.constraintNames{num_constr+1,1} = 'GLCratio';
    NRAmodel.constraintType{num_constr+1,1} = '=';
    NRAmodel.rhs(num_constr+1) = 0;
    NRAmodel.A(num_constr+1,find_cell('RXN_GLCptspp',NRAmodel.varNames)) = 1;
    NRAmodel.A(num_constr+1,find_cell('RXN_GLCtex',NRAmodel.varNames)) = -1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the thermodynamically permissible bounds
% (with or without +-pH allowance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Conc_idx = getVariables(NRAmodel,{'LC'});
Conc_variables = NRAmodel.varNames(Conc_idx);

tempLC_LB = log(1.e-7)*ones(size(Conc_variables,1),1);
tempLC_UB = log(1.e-2)*ones(size(Conc_variables,1),1);
for i=1:size(Conc_variables,1)
    idLC_TFBA = find_cell(Conc_variables(i),options.tFBAmodel.varNames);
    if ~isempty(idLC_TFBA)
        tempLC_LB(i) = options.tFBAmodel.var_lb(idLC_TFBA);
        tempLC_UB(i) = options.tFBAmodel.var_ub(idLC_TFBA);
    end
end

Minmax_ref = [tempLC_LB tempLC_UB];

if isfield(options,'relaxpH') && options.relaxpH~=0
    DGo_idx = getAllVar(NRAmodel,{'DGo'});
    Thermo_rxns = strrep(NRAmodel.varNames(DGo_idx),'DGo_','');
    Thermo_rxns = setdiff(Thermo_rxns,{'Growth';'Htex'});
    DGo_ref = [NRAmodel.var_lb(DGo_idx) NRAmodel.var_ub(DGo_idx)];
    [tempLC_LB,tempLC_UB,DGo_lb,DGo_ub] = getLCvaluesforPH(Conc_variables,options,Minmax_ref,Thermo_rxns,DGo_ref);
    if isfield(options,'DGconstraints') && options.DGconstraints==1
        NRAmodel.var_lb(DGo_idx) = DGo_lb;
        NRAmodel.var_ub(DGo_idx) = DGo_ub;
    end
end

idMets_TFBA = find_cell(strrep(Conc_variables,'LC_',''),options.tFBAmodel.mets);
refLC = options.tFBAmodel.SolLC(idMets_TFBA);
ConcViolations.Conc_variables = Conc_variables;
ConcViolations.bounds = [log([exp(1).^tempLC_LB]./[exp(1).^refLC]) log([exp(1).^tempLC_UB]./[exp(1).^refLC])];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the relaxed bounds for each metabolite
% We allow it to violate up to the physiological bounds defined
% in compartment data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size(Conc_variables,1)
    idLC_TFBA = find_cell(strrep(Conc_variables(i),'LC_',''),options.tFBAmodel.mets);
    if ~isempty(idLC_TFBA)
        LC_comp_symbol = options.tFBAmodel.metCompSymbol(idLC_TFBA);
        tempLC_LB(i) = log(options.tFBAmodel.CompartmentData.compMinConc(find_cell(LC_comp_symbol,options.tFBAmodel.CompartmentData.compSymbolList)));
        tempLC_UB(i) = log(options.tFBAmodel.CompartmentData.compMaxConc(find_cell(LC_comp_symbol,options.tFBAmodel.CompartmentData.compSymbolList)));
    end
end
% Retrieve the ID in NRA model
idLC_NRA = find_cell(Conc_variables,NRAmodel.varNames);
% Apply to concentrations
NRAmodel.var_lb(idLC_NRA) = log([exp(1).^tempLC_LB]./[exp(1).^refLC]);
NRAmodel.var_ub(idLC_NRA) = log([exp(1).^tempLC_UB]./[exp(1).^refLC]);
% Hydrogens need to be relaxed in a different way
H_idx = find_cell(strcat('LC_',{'h_c';'h_p'}),NRAmodel.varNames);
NRAmodel.var_lb(H_idx) = log(1./10);
NRAmodel.var_ub(H_idx) = log(10);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mkdir(NRA_dir)
save([NRA_dir,'NRAmodel.mat'],'NRAmodel')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = 1;
for i = options.violenzymes
    for x = 0:options.violstep:options.max_viol
        
        fprintf('\n --------------------------------------------\n')
        fprintf('Max violations %d, for %i manipulations.',x,i)
        fprintf('\n --------------------------------------------\n\n')
        
        NRAmodel.rhs(find_cell('TOTAL_PERTUNRATIONS',NRAmodel.constraintNames)) = ...
            numel(consnames)-i;
        
        M = NRAmodel;
        
        sol = solveTFAmodelCplex(NRAmodel,options.timelimit);
        
        if isempty(sol.x)
            fprintf('*Checkpoint break for %i manipulations*\n',i)
            break
        end
        
        constraints = [Conc_variables num2cell(ConcViolations.bounds)];
        [~,constraintsToRelax,fc] = debugAdditionOfConstraints_W(M,constraints,x,options.objrxn,options.timelimit,NRA_dir);
        
        ConcViolations.Enzyme_Manipulations(j) = i;
        ConcViolations.Violations(j) = size(constraintsToRelax,1);
        ConcViolations.Concs{j} = constraintsToRelax;
        ConcViolations.Fold_Change(j) = fc;
        
        fprintf('Fold change = %d\n',fc)
        j = j+1;
    end
end

ConcViolations.Enzyme_Manipulations = ConcViolations.Enzyme_Manipulations';
ConcViolations.Violations = ConcViolations.Violations';
ConcViolations.Concs = ConcViolations.Concs';
ConcViolations.Fold_Change = ConcViolations.Fold_Change';
ConcViolations.options = options;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([NRA_dir,'ConcViolations.mat'],'ConcViolations')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
