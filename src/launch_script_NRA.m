% Script for launching the NRA case studies for a given MCA problem

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add relevant repository paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
git_dir =  '/Users/xxx/yyy/'; % path when the mattfa cloned repository is located
cplex_dir =  '/Users/xxx/zzz/'; % path when the CPLEX Studio is located
addpath(genpath(strcat(git_dir,'mattfa/')))
addpath(genpath(strcat(cplex_dir,'CPLEX_Studio1251/'))) % or any other version - adjust

addpath(genpath('./nra_functions/'))

% set a time limit for the solver
options.timelimit = 1800;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Ecoli (or other) case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.modelcase = 'Ecoli';
 
switch options.modelcase
    case 'Ecoli'
        model_name='D1_FDP1';
        options.biomassRxn = {'Ec_biomass_iJO1366_WT_53p95M'};
        options.fixGLCptsppRatio = 1; % extra constrain to fix the flux ratio of glucose to the ref case (~91% - McCloskey data)
end

save_dir = './../modeldata'; % path where modeldata folder is located

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set bounds for variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conc_fold = 2; % this will be overwritten if other option are enabled below (DG constraints)
flux_fold = 10;
enz_fold = 5;

% change cobra solver
changeCobraSolver('cplex_direct','LP')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% option for concatenating the control coefficients and performing PCA with
% respect to flux or conc control coefficients.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.runFirstTime = 0; % Concatenate the FCCs and CCs and compute Mu FCCs
options.runFirstTimeOBJrxn = 0; % Select extremes based on PCA of the CC for new OBJ rxn

options.getMu = 0; % get growth FCCs and consider Mu


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enumerate all alternative solutions for given problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.enumAlternatives = 0;

% Max number of alternative strategies we want to generate
% Note: If number [options.num_enz_perturb] not fixed or if its large, 
% this can take a while as combinatoric explosion is pretty tremendous....
options.noAlt = 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Symmetric bounds options
% (non symmetric accounts for gene knockouts so lower bound is inf-fold)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do we want symmetric enzyme changes?
options.symmEnz = 0;
% Do we want symmetric flux changes?
options.symmFlux = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option for deciding whether to run all CC sets (PCA 19)
% Else we just use the reference CC set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.runAllSaturations = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Automatic setting of model bounds based on options defined above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.metub = log(conc_fold);
options.metlb = log(1/conc_fold); %default symmetric
if options.symmFlux == 1
    options.fluxub = log(flux_fold);
    options.fluxlb = log(1./flux_fold);
elseif options.symmEnz == 0
    options.fluxub = log(flux_fold);
    options.fluxlb = -1000;
end
if options.symmEnz == 1
    options.enzub = log(enz_fold);
    options.enzlb = options.enzub;
elseif options.symmEnz == 0
    options.enzub = log(enz_fold);
    options.enzlb = 1000;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add explicit DG constraints for the concentrations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.DGconstraints = 1;

% Add the preset bounds for concs (compartment data from original model)
options.presetfoldcons = 1; % Note:This overrides the conc bounds

% Option for adding lb and ub from the tfa model
% Add them if you want to constrain the space a priori
options.boundsTFA = 0; % Note:This overrides the conc bounds

if options.presetfoldcons == 1 && options.boundsTFA == 1
    disp('options.presetfoldcons and options.boundsTFA cannot both be active! Check!')
    keyboard
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adds the additional constraint of x-fold to the conc bounds.
% Put to zero to have it inactive
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.maxFoldConc = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting of OTHER NON-OPTIMIZATION objectives we might want to constrain
% Simple for manual setting of Fluxes/Concs bounds as vectors
% Experimental for reading data from mat files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting of objective options
switch options.modelcase
    case 'Ecoli'
        options.objvars = {'RXN_GLCptspp'}; % define rxn names as 'RXN_RXNNAME'
        options.objvlb = options.fluxlb; % define bounds
        options.objvub = options.fluxub; % define bounds
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set which CC column we want to use for PCA
% usually the same as objective
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.objrxnPCA = {'RXN_GLCptspp'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set fixed number of enzymatic manipulations
% Put zero for no constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.num_enz_perturb = 5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set objective type (Rxn or Enzymes)
% 'Rxn' to use the objrxn & objtype as objectives
% 'Enzymes' minimizes the no of manipulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.setobjective = 'Rxn';

if strcmp(options.setobjective,'Rxn')
    options.objrxn = options.objrxnPCA; % by default use the same rxn for PCA as the objective
    options.objtype = -1;
elseif strcmp(options.setobjective,'Enzymes')
    options.num_enz_perturb = 0; % we don't want to fix the number of manipulations
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concentration violation study option
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.ConcViolationStudy = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directory where to store the data output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NRA_dir = ['./dataNRA_enz',num2str(exp(options.enzub)),...
    'conc',num2str(exp(options.metub)),'flux',num2str(...
    exp(options.fluxub)),'_',char(strrep(options.objrxn,'RXN_','')),'/'];
mkdir(NRA_dir);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Automatic loading of model data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([save_dir,'metSampData',model_name,'.mat'],'tFBAmodel')
options.tFBAmodel = tFBAmodel;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     RUN     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the 19 extreme CC sets from the full CC set.This will save the CCs
% that were chosen based on PCA so we don't have all CCs locally.
if options.runFirstTime == 1
    concatenateCCs_PCA_Stats(save_dir,'Cpagg',model_name)
    concatenateCCs_PCA_Stats(save_dir,'Cxagg',model_name)
    if options.getMu == 1
        getCCforMu(save_dir,model_name,options)
    end
end

if options.runFirstTimeOBJrxn == 1
    fetchCCfromPCA(save_dir,model_name,options);
end

if options.ConcViolationStudy
    options.violenzymes = [1,2,4,7]; % DEFINE NUMBER OF ENZYMES
    options.max_viol = 20; % Max violations (starts from 0)
    options.violstep = 2; % step for the iterations
    % Option to include a positive or negative pH relaxation in the
    % study - it will recalculate DGs and give flexibility to concs
    % NOTE: Some parts for this option are hardcoded for E.coli
    options.relaxpH = 0;
    options.pHdiff = 0.5; % by how much to relax the pH
    options.pHdirection = '+'; % '+' or '-' directions
    load([save_dir,model_name,'.mat'])
    eval(['options.ThermoModel = ',model_name,';']);
    case_dir = [NRA_dir,'ConcViolationsStudy/'];
    ConcViolations = runConcViolationStudies(save_dir,model_name,options,case_dir);
else
    % Normal NRA run
    run_NRA(save_dir,model_name,options,NRA_dir)
end






