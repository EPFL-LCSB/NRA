function [NRAmodel, ind_bfuse] = prepMCAcoefficientsForNRA(controlflux, controlconc, consnames,kinModel, options, fluxnames )
% function that builds the NRA model

if nargin<6
    fluxnames = consnames;
end

% Build the N matrix refered to as K in the NRA
N = [-eye(size(controlconc,1)) zeros(size(controlconc,1),size(controlflux,1)) controlconc -controlconc;...
    zeros(size(controlflux,1),size(controlconc,1)) -eye(size(controlflux,1)) controlflux -controlflux];
% Variable names of u
v = [strcat('LC_',kinModel.metabolites); strcat('RXN_',fluxnames); strcat('F_',consnames); strcat('R_', consnames)];
% RHS of zeros
rhs = zeros(size(N,1),1);
constype(1:size(N,1),1) = {'='};
vartype1(1:size(controlconc,1),1) = {'C'};
vartype2(1:size(controlflux,1),1) = {'C'};
vartype3(1:size(consnames),1) = {'C'};
vartype4(1:size(consnames),1) = {'C'};
vartype = [vartype1; vartype2; vartype3; vartype4];
model.A = N;
model.constraintNames = [kinModel.metabolites; fluxnames];
model.varNames = v;
model.constraintType = constype;
model.vartypes = vartype;
model.rhs = rhs;
model.description = 'control';

% Fixing the bounds for fluxes and concentrations
model.var_ub = [options.metub*ones(length(kinModel.metabolites),1); options.fluxub*ones(length(fluxnames),1);...
    options.enzub*ones(length(consnames),1); options.enzlb*ones(length(consnames),1)];
model.var_lb = [options.metlb*ones(length(kinModel.metabolites),1); options.fluxlb*ones(length(fluxnames),1);...
    zeros(length(consnames),1); zeros(length(consnames),1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BINARY CONTROL ON ENZYMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[~,num_vars] = size(model.A);
forward=getAllVar(model, {'F'});
backward=getAllVar(model, {'R'});
% Add binary variable BFUSE for whether or not we modify an enzyme
for i=1:length(forward)
    model.varNames{length(model.varNames)+1} = strcat('BFUSE_', consnames{i});
    model.var_ub(length(model.varNames)) = 1;
    model.var_lb(length(model.varNames)) = 0;
    model.vartypes(length(model.varNames)) = {'B'};
    model.f(length(model.varNames)) = 1;
end
% Add constraint BFMER_:  F_ + B_+ 10000*BFUSE < 10000
for i=1:length(forward)
    [num_constr,~] = size(model.A);
    model.rhs(num_constr+1,1) = 10000;
    model.constraintNames{num_constr+1,1} = strcat('BFMER_',num2str(i));
    model.constraintType{num_constr+1,1} = '<';
    model.A(num_constr+1,forward(i)) = 1;
    model.A(num_constr+1,backward(i)) = 1;
    model.A(num_constr+1,num_vars+i) = 10000;
    ind_bfuse(i,1) = num_vars+i;
end
% Add the binaries FU_ and BU_ to allow up or down regulation of an enzyme
for i=1:length(forward)
    model.varNames{length(model.varNames)+1,1} = strcat('FU_',consnames{i});
    model.var_ub(length(model.varNames),1) = 1;
    model.var_lb(length(model.varNames),1) = 0;
    model.vartypes{length(model.varNames),1} = 'B';
    model.varNames{length(model.varNames)+1,1} = strcat('BU_',consnames{i});
    model.var_ub(length(model.varNames),1) = 1;
    model.var_lb(length(model.varNames),1) = 0;
    model.vartypes{length(model.varNames),1} = 'B';
end

% Find FU_ and BU_ indices so that we can add contraints
% such that an enzyme cannot be both up and down regulated
% at the same time
[~,indfu] = ismember(strcat('FU_',consnames), model.varNames);
[~,indbu] = ismember(strcat('BU_',consnames), model.varNames);
% We then add the constraint SU_: FU_ + BU_ < 1
for i=1:length(forward)
    [num_constr,~] = size(model.A);
    model.rhs(num_constr+1,1) =  1;
    model.constraintNames{num_constr+1,1} = strcat('SU_',consnames{i});
    model.constraintType{num_constr+1,1} = '<';
    model.A(num_constr+1,indfu(i)) = 1;
    model.A(num_constr+1,indbu(i)) = 1;
end

% Find F_ and R_
[~,indf] = ismember(strcat('F_',consnames), model.varNames);
[~,indr] = ismember(strcat('R_',consnames), model.varNames);
% Add constraints UF_ and UB_
% UF_: F_ -10000 FU_ < 0
% UB_: B_ -10000 BU_ < 0
for i=1:length(indf)
    [num_constr,~] = size(model.A);
    model.rhs(num_constr+1,1) =  0;
    model.constraintNames{num_constr+1,1} = strcat('UF_',consnames{i});
    model.constraintType{num_constr+1,1} = '<';
    model.A(num_constr+1,indf(i)) = 1;
    model.A(num_constr+1,indfu(i)) = -10000;
    [num_constr,~] = size(model.A);
    model.rhs(num_constr+1,1) =  0;
    model.constraintNames{num_constr+1,1} = strcat('UB_',consnames{i});
    model.constraintType{num_constr+1,1} = '<';
    model.A(num_constr+1,indr(i)) = 1;
    model.A(num_constr+1,indbu(i)) = -10000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DG CONSTRAINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(options,'DGconstraints') && options.DGconstraints==1
    Tmodel = options.tFBAmodel;
    Thermo_rxns = strrep(Tmodel.constraintNames(getAllCons(Tmodel,{'G'})),'G_','');
    
    % Add extracellular metabolite concentrations
    ext_mets = Tmodel.mets(find_cell('e',Tmodel.metCompSymbol));
    ext_mets = setdiff(ext_mets,{'h2o_e','h_e'});
    for j=1:size(ext_mets,1)
        [num_constr,num_vars] = size(model.A);
        model.varNames{num_vars+1,1} = strcat('LC_', ext_mets{j});
        model.var_lb(num_vars+1,1)   = options.metlb;
        model.var_ub(num_vars+1,1)   = options.metub;
        model.vartypes{num_vars+1,1} = 'C';
        model.A(num_constr,num_vars+1) = 0;
    end
    
    for i = 1:size(Thermo_rxns,1);
        if strcmp(Thermo_rxns{i},'Growth') || strcmp(Thermo_rxns{i},'Htex')
            continue
        end
        % Grab the reactants & coeffs from tfa model
        rxn_direction = sign(Tmodel.SolFluxes(find_cell(Thermo_rxns(i),Tmodel.rxns)));
        index = find_cell(strcat('G_',Thermo_rxns(i)),Tmodel.constraintNames);
        variables = Tmodel.varNames(find(Tmodel.A(index,:)));
        coefficients = full(Tmodel.A(index,find(Tmodel.A(index,:))));
        
        % Keep only the concentrations from the model DG constraints
        b = find(~cellfun(@isempty,strfind(variables,'LC_')));
        variables = variables(b);
        coefficients = coefficients(b);
        
        % Remove H2O from the reactants
        c = find_cell('cpd00001',Tmodel.metSEEDID(find_cell(strrep(variables,'LC_',''),Tmodel.mets)));
        variables(c) = [];
        coefficients(c) = [];
        
        % Get DGo bounds
        c = find_cell(strcat('DGo_', Thermo_rxns{i}),Tmodel.varNames);
        if isempty(c)
            c = find_cell(strcat('DGMo_', Thermo_rxns{i}),Tmodel.varNames);
        end
        DGo_lb = Tmodel.var_lb(c);
        DGo_ub = Tmodel.var_ub(c);
        
        % Add DGo variables
        [num_constr,num_vars] = size(model.A);
        model.varNames{num_vars+1,1} = strcat('DGo_', Thermo_rxns{i});
        model.var_lb(num_vars+1,1)   = DGo_lb;
        model.var_ub(num_vars+1,1)   = DGo_ub;
        model.vartypes{num_vars+1,1} = 'C';
        model.A(num_constr,num_vars+1) = 0;
        
        DGo_index = size(model.varNames,1);
        LC_indexes = find_cell(variables,model.varNames);
        
        % DG: DGo + RT*StoichCoefProd1*M_prod1 + RT*StoichCoefSub1*M_subs1
        %          ...< - (RT*StoichCoefSub1*log(ref_subs1) + RT*StoichCoefSub2*log(ref_subs2))'
        % Formulate the constraint
        CLHS.varIDs    = [DGo_index  LC_indexes  ];
        CLHS.varCoeffs = [1          coefficients];
        
        RHS = 0;
        for k = 1:size(variables,1)
            var_ref = Tmodel.SolLC(find_cell(strrep(variables(k),'LC_',''),Tmodel.mets));
            RHS = RHS + coefficients(k)*var_ref;
        end
        
        [num_constr,~] = size(model.A);
        model.constraintNames{num_constr+1,1} = strcat('DG_', Thermo_rxns{i});
        if rxn_direction>0
            model.constraintType{num_constr+1,1}  = '<';
        else
            model.constraintType{num_constr+1,1}  = '>';
        end
        model.rhs(num_constr+1) = -RHS;
        model.A(num_constr+1,CLHS.varIDs) = CLHS.varCoeffs;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIAL CONCENTRATION BOUNDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(options,'boundsTFA') && options.boundsTFA==1
    if isfield(options,'DGconstraints') && options.DGconstraints==1
        tmp_mets = setdiff(options.tFBAmodel.mets,{'h2o_c';'h2o_p';'h2o_e';'h_e'});
    else
        tmp_mets = kinModel.metabolites;
    end
    tempLC_LB = log(1.e-7)*ones(size(tmp_mets,1),1);
    tempLC_UB = log(1.e-2)*ones(size(tmp_mets,1),1);
    for i=1:size(tmp_mets,1)
        idLC_TFBA = find_cell(strcat('LC_',tmp_mets(i)),options.tFBAmodel.varNames);
        if ~isempty(idLC_TFBA)
            tempLC_LB(i) = options.tFBAmodel.var_lb(idLC_TFBA);
            tempLC_UB(i) = options.tFBAmodel.var_ub(idLC_TFBA);
        end
    end
    idMets_TFBA = find_cell(tmp_mets,options.tFBAmodel.mets);
    refLC = options.tFBAmodel.SolLC(idMets_TFBA);
    % Retrieve the ID in NRA model
    idLC_NRA = find_cell(strcat('LC_',tmp_mets),model.varNames);
    % Apply to concentrations
    model.var_lb(idLC_NRA) = log([exp(1).^tempLC_LB]./[exp(1).^refLC]);
    model.var_ub(idLC_NRA) = log([exp(1).^tempLC_UB]./[exp(1).^refLC]);
    if isfield(options,'maxFoldConc') && options.maxFoldConc~=0
        % Update LB
        my_tempLB = [exp(1).^tempLC_LB]./[exp(1).^refLC];
        isLess = my_tempLB<1/options.maxFoldConc;
        my_tempLB(isLess) = 1/options.maxFoldConc;
        % Update UB
        my_tempUB = [exp(1).^tempLC_UB]./[exp(1).^refLC];
        isMore = my_tempUB>options.maxFoldConc;
        my_tempUB(isMore) = options.maxFoldConc;
        % update bounds so our TFA does not exceed desired physiology
        model.var_lb(idLC_NRA) = log(my_tempLB);
        model.var_ub(idLC_NRA) = log(my_tempUB);
    end
elseif isfield(options,'presetfoldcons') && options.presetfoldcons==1
    if isfield(options,'DGconstraints') && options.DGconstraints==1
        tmp_mets = setdiff(options.tFBAmodel.mets,{'h2o_c';'h2o_p';'h2o_e';'h_e'});
    else
        tmp_mets = kinModel.metabolites;
    end
    tempLC_LB = log(1.e-7)*ones(size(tmp_mets,1),1);
    tempLC_UB = log(1.e-2)*ones(size(tmp_mets,1),1);
    for i=1:size(tmp_mets,1)
        idLC_TFBA = find_cell(tmp_mets(i),options.tFBAmodel.mets);
        if ~isempty(idLC_TFBA)
            LC_comp_symbol = options.tFBAmodel.metCompSymbol(idLC_TFBA);
            tempLC_LB(i) = log(options.tFBAmodel.CompartmentData.compMinConc(find_cell(LC_comp_symbol,options.tFBAmodel.CompartmentData.compSymbolList)));
            tempLC_UB(i) = log(options.tFBAmodel.CompartmentData.compMaxConc(find_cell(LC_comp_symbol,options.tFBAmodel.CompartmentData.compSymbolList)));
        end
    end
    idMets_TFBA = find_cell(tmp_mets,options.tFBAmodel.mets);
    refLC = options.tFBAmodel.SolLC(idMets_TFBA);
    % Retrieve the ID in NRA model
    idLC_NRA = find_cell(strcat('LC_',tmp_mets),model.varNames);
    % Apply to concentrations
    model.var_lb(idLC_NRA) = log([exp(1).^tempLC_LB]./[exp(1).^refLC]);
    model.var_ub(idLC_NRA) = log([exp(1).^tempLC_UB]./[exp(1).^refLC]);
    if isfield(options,'maxFoldConc') && options.maxFoldConc~=0
        % Update LB
        my_tempLB = [exp(1).^tempLC_LB]./[exp(1).^refLC];
        isLess = my_tempLB<1/options.maxFoldConc;
        my_tempLB(isLess) = 1/options.maxFoldConc;
        % Update UB
        my_tempUB = [exp(1).^tempLC_UB]./[exp(1).^refLC];
        isMore = my_tempUB>options.maxFoldConc;
        my_tempUB(isMore) = options.maxFoldConc;
        % update bounds so our TFA does not exceed desired physiology
        model.var_lb(idLC_NRA) = log(my_tempLB);
        model.var_ub(idLC_NRA) = log(my_tempUB);
    end
    % Hydrogens cannot change their concentration IF we account for
    % thermodynamics (else we set if to the conc bound)
    H_idx = find_cell(strcat('LC_',{'h_c';'h_p'}),model.varNames);
    if isfield(options,'DGconstraints') && options.DGconstraints==1
        model.var_lb(H_idx) = 0;
        model.var_ub(H_idx) = 0;
    else
        model.var_lb(H_idx) = options.metlb;
        model.var_ub(H_idx) = options.metub;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BINARY CONTROL ON ENZYMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model.f = zeros(length(model.varNames),1);
model.f(ind_bfuse) = 1;
% Here it is maximization but actually it minimizes the total use of the
% sum of enzymes that change
model.objtype = -1;

% Constrain to fix number of enzymatic pertuNRAtions
if isfield(options,'num_enz_perturb') && options.num_enz_perturb~=0
    [num_cons,num_vars] = size(model.A);
    new_constraint = zeros(1,num_vars);
    new_constraint(ind_bfuse) = 1;
    model.A(num_cons+1,:) = new_constraint;
    model.rhs(num_cons+1) = numel(consnames)-options.num_enz_perturb;
    model.constraintNames{num_cons+1} = ['TOTAL_PERTUNRATIONS'];
    model.constraintType{num_cons+1} = '=';
end

NRAmodel = model;
end