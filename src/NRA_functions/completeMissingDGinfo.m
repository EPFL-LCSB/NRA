function DGCompleteModel = completeMissingDGinfo(model, verboseFlag, flagToAddLnThermoDisp)

if ~exist('verboseFlag','var')
    verboseFlag = 1;
end

if ~exist('flagToAddLnThermoDisp','var') || isempty(flagToAddLnThermoDisp)
    flagToAddLnThermoDisp = false;
end

if strcmp(model.thermo_units,'kJ/mol')
    GAS_CONSTANT = 8.314472/1000; % kJ/(K mol)
else
    GAS_CONSTANT = 1.9858775/1000; % Kcal/(K mol)
    % value for the bigM in THermo constraints. This will also be the bound value
    bigMtherm = 1.e+6;
    DGR_lb = -bigMtherm; %kcal/mol
    DGR_ub =  bigMtherm; %kcal/mol
end
TEMPERATURE = 298.15; % K
RT = GAS_CONSTANT*TEMPERATURE;

% value for the bigM in big M constraints such as:
% UF_rxn: F_rxn - M*FU_rxn < 0
bigM = 1000;

% ---- FOR EACH METABOLITE (THAT HAS NO THERMODYNAMIC INFORMATION) ----
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ->->->->->-> METABOLITES WITHOUT LC-variables ->->->->->->->
% Ids of LC(Log-concentration) variables in the model
LC_id = getAllVar(model,{'LC'});
% Metabolite names of metabolites that have a corresponding LC_variable in the (tFBA-)model
LCmetNames = regexprep(model.varNames(LC_id),'^(LC_)','');
% Metabolite names of metabolites that do NOT have a corresponding LC_variable in the (tFBA-)model
NoLCmetNames = setdiff(model.mets,LCmetNames);
% Get the indices of these metabolites with missing LC
idNoLCmetNames = find_cell(NoLCmetNames,model.mets);
if size(idNoLCmetNames,2) ~= size(NoLCmetNames,1)
    error('Inconsistent number of missing LC variables')
end
% <-<-<-<-<-<- METABOLITES WITHOUT LC-variables <-<-<-<-<-<-<-
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


for i = 1:size(idNoLCmetNames,2)
    idMetIntFBA = idNoLCmetNames(i);
    % exclude protons and water and those without deltaGF
    % P_met: P_met - RT*LC_met + DGFerr_met = DGF_met
    metDeltaGF = model.metDeltaGFtr(idMetIntFBA);
    metComp = model.metCompSymbol{idMetIntFBA};
    Comp_index = find(ismember(model.CompartmentData.compSymbolList,metComp));
    metLConc_lb = log(model.CompartmentData.compMinConc(Comp_index));
    metLConc_ub = log(model.CompartmentData.compMaxConc(Comp_index));
    if verboseFlag
        fprintf('generating thermo variables for %s\n',model.mets{idMetIntFBA});
    end
    if  (metDeltaGF > 1E6)
        model = addNewVariableInTFA(model,strcat('LC_',model.mets{idMetIntFBA}),'C',[metLConc_lb metLConc_ub]);
    else
        error('How did this metabolite end up without DG?')
    end
end

% ---- FOR EACH REACTION (THAT HAS NO THERMODYNAMIC INFORMATION) ----
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ->->->->->-> REACTIONS WITHOUT THERMO ->->->->->->->
% The reactions that have NO thermodynamic information are:
idRxnsWithoutThermo = find(model.rxnThermo==0);
% From these reactions we should exclude
% (1) The biomass reaction, because it is not stoichiometrically balanced,
% so the thermodynamic constraints would be just wrong:
idBiomass = find_cell(cellfun(@(x) (x(3:end)),model.varNames(model.f==1),'UniformOutput',false),model.rxns);
idRxnsWithoutThermo = setdiff(idRxnsWithoutThermo,idBiomass);
% (2) The boundary fluxes (so-called drain fluxes in the LCSB-lingo)
idDrains = find(sum(model.S~=0,1) == 1);
idRxnsWithoutThermo = setdiff(idRxnsWithoutThermo,idDrains);
% Keep a backup from the ractions that had no thermo
model.rxnThermoBackup = model.rxnThermo;
% <-<-<-<-<-<- REACTIONS WITHOUT THERMO <-<-<-<-<-<-<-
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

disp('Generating thermodynamic constraints for reactions');
for i = 1:size(idRxnsWithoutThermo,1)
    idRxnIntFBA = idRxnsWithoutThermo(i);
    % Get the net flux index of this reaction
    NF_flux_index = find(ismember(model.varNames,strcat('NF_',model.rxns{idRxnIntFBA})));
    % Get the forward and reverse flux index of this reaction
    F_flux_index = find(ismember(model.varNames,strcat('F_',model.rxns{idRxnIntFBA})));
    R_flux_index = find(ismember(model.varNames,strcat('R_',model.rxns{idRxnIntFBA})));
    % Extract the reaction directionality in this particular FDP
    if model.var_lb(NF_flux_index)>0 && model.var_ub(NF_flux_index)>0
        CurDirection = 1;
    elseif model.var_lb(NF_flux_index)<0 && model.var_ub(NF_flux_index)<0
        CurDirection = -1;
    else
        error('something is wrong with the current FDP')
    end
    % We assign the DGerror of the reaction as 2 kcal/mol, as an
    % approximation to the median error assigned to the rest reactions
    % median(model.rxnDeltaGRerr(model.rxnDeltaGRerr~=10000000))
    % ans = 2
    RxnDGerror  = 2;
    model.rxnDeltaGRerr(idRxnIntFBA) = RxnDGerror;
    % Similarly we assign the DG naught of the reaction as +/- 2.5 kcal/mol
    % median(model.rxnDeltaGR(model.rxnDeltaGR~=10000000))
    % ans= -2.5927
    %     RxnDGnaught = -2.5*CurDirection;
    % We assign a high value for now, but we'll have to define it more
    % systematically later on:
    RxnDGnaught = -100*CurDirection;
    model.rxnDeltaGR(idRxnIntFBA) = RxnDGnaught;
    
    if verboseFlag
        fprintf('generating thermo constraint for %s\n',model.rxns{idRxnIntFBA});
    end
    
    model = addNewVariableInTFA(model,strcat('DGM_',model.rxns{idRxnIntFBA}),'C',[DGR_lb DGR_ub]);
    DGM_index = size(model.varNames,1);
    model = addNewVariableInTFA(model,strcat('DGMo_',model.rxns{idRxnIntFBA}),'C', RxnDGnaught + [-RxnDGerror RxnDGerror]);
    DGMo_index = size(model.varNames,1);
    
    % GF Dec 2015: Create DG constraint:
    % G: -DGM_rxn + DGMo_Rxn + RT*LC_prod1 + RT*LC_prod2 + ... - RT*LC_subs1 - RT*LC_subs2 - ... = - DGMo_Rxn'
    % Find substrates and products of the reaction
    id_prods = model.S(:,idRxnIntFBA)>0;
    id_subs  = model.S(:,idRxnIntFBA)<0;
    RxnProducts   = model.mets(id_prods);
    RxnSubstrates = model.mets(id_subs);
    % Get the corresponding LC_-variables
    [LogicalId_ProdsLCs_InVarNames,Id_ProdsLCs_InVarNames] = ismember(strcat('LC_',RxnProducts),model.varNames);
    LC_Prods_indexes = Id_ProdsLCs_InVarNames(LogicalId_ProdsLCs_InVarNames);
    [LogicalId_SubsLCs_InVarNames,Id_SubsLCs_InVarNames] = ismember(strcat('LC_',RxnSubstrates),model.varNames);
    LC_Subs_indexes  = Id_SubsLCs_InVarNames(LogicalId_SubsLCs_InVarNames);
    % Formulate the constraint
    CLHS.varIDs    = [DGM_index  DGMo_index LC_Prods_indexes'                    LC_Subs_indexes'                    ];
    CLHS.varCoeffs = [-1         1          RT*ones(1,size(LC_Prods_indexes,1))  -RT*ones(1,size(LC_Subs_indexes,1)) ];
    model = addNewConstraintInTFA(model,strcat('G_',model.rxns{idRxnIntFBA}),'=',CLHS, 0);
    
    % create the use variables constraints and connect them to the
    % deltaG if the reaction has thermo constraints
    % FU_rxn: 1000 FU_rxn + DGM_rxn < 1000 - epsilon
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Note: The reason why we need this epsilon is because when the    %
    % FU_rxn binary variable is equal to one (i.e. FU_rxn = 1) then    %
    % what is left is that DGR_rxn < 0. However, CPLEX always uses <=, %
    % meaneng that CPLEX can find a feasible solution for DGR_rxn = 0. %
    % Therefore we use epsilon as something very small to assure that  %
    % the constraint is fullfilling its purpose. Epsilon is used for   %
    % exactly the same reason in the symmetric constraint later on.    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    epsilon = 1e-6;
    
    % Check if FU variable exists. If not, then create it.
    FUvarName = strcat('FU_',model.rxns{idRxnIntFBA});
    FU_index = find_cell(FUvarName,model.varNames);
    if isempty(FU_index)
        model = addNewVariableInTFA(model,FUvarName,'B',[0 1]);
        FU_index = size(model.varNames,1);
    end
    %     if (model.rxnThermo(idRxnIntFBA) == 1) % WE HAVE TO SKIP THIS
    CLHS.varIDs    = [DGM_index FU_index];
    CLHS.varCoeffs = [1        +bigMtherm];
    model = addNewConstraintInTFA(model,FUvarName,'<',CLHS, bigMtherm - epsilon);
    %     end
    % BU_rxn: 1000 BU_rxn - DGR_rxn < 1000 + epsilon
    % Check if BU variable exists. If not, then create it.
    BUvarName = strcat('BU_',model.rxns{idRxnIntFBA});
    BU_index = find_cell(BUvarName,model.varNames);
    if isempty(BU_index)
        model = addNewVariableInTFA(model,BUvarName,'B',[0 1]);
        BU_index = size(model.varNames,1);
    end
    %     if (model.rxnThermo(idRxnIntFBA) == 1)
    CLHS.varIDs    = [DGM_index BU_index];
    CLHS.varCoeffs = [-1       +bigMtherm];
    model = addNewConstraintInTFA(model,BUvarName,'<',CLHS, bigMtherm + epsilon);
    %     end
    % create the prevent simultaneous use constraints
    % U_rxn: FU_rxn + BU_rxn <= 1
    SUconstName = strcat('SU_',model.rxns{idRxnIntFBA});
    if isempty(find_cell(SUconstName,model.constraintNames))
        CLHS.varIDs    = [FU_index BU_index];
        CLHS.varCoeffs = [+1       +1];
        model = addNewConstraintInTFA(model,SUconstName,'<',CLHS, 1);
    end
    % create constraints that control fluxes with their use variables
    % UF_rxn: F_rxn - 100 FU_rxn < 0
    UFconstName = strcat('UF_',model.rxns{idRxnIntFBA});
    if isempty(find_cell(UFconstName,model.constraintNames))
        CLHS.varIDs    = [F_flux_index FU_index];
        CLHS.varCoeffs = [+1           -bigM];
        model = addNewConstraintInTFA(model,UFconstName,'<',CLHS, 0);
    end
    % UR_rxn: R_rxn - 100 RU_rxn < 0
    URconstName = strcat('UR_',model.rxns{idRxnIntFBA});
    if isempty(find_cell(URconstName,model.constraintNames))
        CLHS.varIDs    = [R_flux_index BU_index];
        CLHS.varCoeffs = [+1           -bigM];
        model = addNewConstraintInTFA(model,URconstName,'<',CLHS, 0);
    end
    
    
    % Below we add to the model a variable for the log of the thermodynamic displacement (Gamma)
    %
    % Gamma = (1/Keq)*(C_Prod1^StoichCoeffProd1 * ...)/(C_Sub1^StoichCoeffSub1 * ... )
    % Keq   = exp(-DGnaught/RT)
    %
    % Combining the two equations above, formulate the constraint:
    % LnGamma: ln(Gamma) + StoichCoefSubs1 * LCsubs1 + StoichCoefSubs2 * LCsubs2 + ...
    %                    - StoichCoefProd1 * LCprod1 - StoichCoefProd2 * LCprod2 - ...
    %                    - (1/RT)*DGo_Rxn = 0
    % NOTE1: The formulation above should be exactly equivalent to:
    % ln(Gamma) - (1/RT)*DG_rxn = 0
    % BUT, the reason why the formulation above is not 100% equivalent
    % to the initial formulation is due to the DG of the transport
    % fluxes that is not added explicitly, but only through the DG.
    % Therefore We adopt the latter formulation, that is also simpler.
    
    if flagToAddLnThermoDisp
        LNvarName = strcat('LnGamma_',model.rxns{idRxnIntFBA});
        if isempty(find_cell(LNvarName,model.varNames))
            model = addNewVariableInTFA(model,strcat('LnGamma_',model.rxns{idRxnIntFBA}),'C',[-1000 1000]);
            LnGamma_index = size(model.varNames,1);
            
            CLHS.varIDs    = [LnGamma_index     DGM_index];
            CLHS.varCoeffs = [1                 -1/RT ];
            model = addNewConstraintInTFA(model, strcat('ThermoDisp_', model.rxns{idRxnIntFBA}), '=', CLHS, 0);
            
            % Here we add MCA constraints on LnGamma_ so that we ensure a
            % rection is not at equilibrium ie. Gamma ~= 1 !!!
            % We add the following constraints:
            % FU_ThermoDisp_: ln(Gamma) + 1000*FU_ <  1000 + Epsilon1
            % BU_ThermoDisp_: ln(Gamma) - 1000*BU_ > -1000 + Epsilon2
            % We want to bound Gamma to be below 0.99 for the forward
            % reaction and it should be above corresponding 1.0101 for
            % the backward reaction (1/0.99=1.0101). REmember we are
            % binding the log of gamma so we take the log of the
            % thermodynamic displacements.
            Epsilon1 = log(0.99);
            Epsilon2 = log(1.0101);
            
            CLHS.varIDs    = [LnGamma_index     FU_index];
            CLHS.varCoeffs = [1                 1000    ];
            model = addNewConstraintInTFA(model, strcat('FU_ThermoDisp_', model.rxns{idRxnIntFBA}), '<', CLHS, 1000 + Epsilon1);
            
            CLHS.varIDs    = [LnGamma_index     BU_index];
            CLHS.varCoeffs = [1                 -1000   ];
            model = addNewConstraintInTFA(model, strcat('BU_ThermoDisp_', model.rxns{idRxnIntFBA}), '>', CLHS, -1000 + Epsilon2);
        end
    end
    % We added the thermodynamic information for this reaction
    model.rxnThermo(idRxnIntFBA) = 1;
    
    % Adjusting the objective
    IDobjective = find(model.f==1);
    model.f = zeros(size(model.A,2),1);
    model.f(IDobjective) = 1;
    
end

% In the next step we are going to sample the concentrations. Based on
% these concentrations we will obtain consistent reaction-DG values and
% reaction-DGerror values. Here we have to assure that these reaction-DG
% values are consistent with the directionality of the current FDP. To
% check this we can run a min-max on the newly added reaction-DG values,
% and make sure that they agree with the FDP for their entire range. If
% not, we would have to assign higher absolut value of DG-naught in the
% desired direction
Check_DG_Flux_Consistency(model)

% Keep the indices of the reactions got thermo information in this
% function:
model.rxnsHadNoThermo = zeros(size(model.rxns));
model.rxnsHadNoThermo(idRxnsWithoutThermo) = 1;

% Set the model to the output
DGCompleteModel = model;

end

function Check_DG_Flux_Consistency(model)
% Here we have to assure that these reaction-DG
% values are consistent with the directionality of the current FDP. To
% check this we can run a min-max on the newly added reaction-DG values,
% and make sure that they agree with the FDP for their entire range. If
% not, we would have to assign higher absolut value of DG-naught in the
% desired direction

% Get all the indices of the newly added reaction-DG values ('DGM_*')
DGMids = getAllVar(model,{'DGM'});
% Run thermominmax for these variables
mmDGMids = runTMinMax(model,model.varNames(DGMids));
% Find the indices of the corresponding NF-variables for these reactions
idCorrespNF = find_cell(regexprep(model.varNames(DGMids),'^(DGM_)','NF_'),model.varNames);
% Get the range of these NF-variables in the model
lbub_NF = [model.var_lb(idCorrespNF) model.var_lb(idCorrespNF)];
% Check consistency:
for i = 1:size(mmDGMids,1)
    if lbub_NF(i,1)>0 && lbub_NF(i,2)>0
        if ~ (mmDGMids(i,1)<0 && mmDGMids(i,2)<0)
            fprintf('NF bounds: %s %s\n',num2str(lbub_NF(i,1)),num2str(lbub_NF(i,2)))
            fprintf('DG bounds: %s %s\n',num2str(mmDGMids(i,1)),num2str(mmDGMids(i,2)))
            error(['Reaction DG of ',model.varNames{DGMids(i)},' is inconsistent!'])
        end
    elseif lbub_NF(i,1)<0 && lbub_NF(i,2)<0
        if ~ (mmDGMids(i,1)>0 && mmDGMids(i,2)>0)
            fprintf('NF bounds: %s %s\n',num2str(lbub_NF(i,1)),num2str(lbub_NF(i,2)))
            fprintf('DG bounds: %s %s\n',num2str(mmDGMids(i,1)),num2str(mmDGMids(i,2)))
            error(['Reaction DG of ',model.varNames{DGMids(i)},' is inconsistent!'])
        end
    else
        error('We have not locked the model to an FDP')
    end
end

end