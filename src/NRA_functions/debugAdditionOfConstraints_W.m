function [constraintsToApply,constraintsToRelax,fc] = debugAdditionOfConstraints_W(model,ConstraintTable,violations,objrxn,timelimit,NRA_dir)
% When introducing constraints, find which ones have to be relaxed in order
% to get a feasible solution.
%
% INPUTS:
% - model: Input TFA model
% - ConstraintTable: Set of variable names for Log concentrations, or
%         reaction fluxes. This should be a cell with three columns:
%         1st column: variable names
%         2nd column: lower bound
%         3rd column: upper bound
%
% OUTPUTS:
% - constraintsToApply: Which constraints can be immediately applied
% - constraintsToRelax: which of the constraints should be relaxed

% Initialize a variable to store all the indices of the new variables for
% debugging
IndOfDebugVars = [];
% model.var_lb(find(model.f)) = 4.1; %growth rate, you can test different growth rates
% Reset the vector of the objective function (zeros)
model.f = zeros(length(model.f),1);
% Find the indices of the variables to be constrained in the varNames
[~, Constr_idx] = ismember(ConstraintTable(:,1), model.varNames);
[~,num_vars] = size(model.A);

% Some contraints of the conc not present in smaller models so remove zeros
Constr_idx = Constr_idx(find(Constr_idx));
[~, id_conc_add] = ismember(model.varNames(Constr_idx),ConstraintTable(:,1));
ConstraintTable = ConstraintTable(id_conc_add,:);

% Create new binary variable DebugRelax,
for i=1:length(ConstraintTable(:,1))
    model.varNames{num_vars+i,1} = strcat('DebugRelax_',model.varNames{Constr_idx(i)});
    model.var_ub(num_vars+i,1) = 1;
    model.var_lb(num_vars+i,1) = 0;
    model.vartypes{num_vars+i,1} = 'B';
    % set the vector of the objective function for the newly generated
    % variables to 1
    model.f(num_vars+i,1) = 1;
    IndOfDebugVars(i,1) = num_vars+i;
end

% Create a new constraint for the newly generated variable
[num_constr,~] = size(model.A);
for i=1:length(ConstraintTable(:,1))
    model.constraintNames{num_constr+i} = strcat('DebugRelax_ub_', model.varNames{Constr_idx(i)});
    model.rhs(num_constr+i) = roundsd(ConstraintTable{i,3},5);
    model.A(num_constr+i, Constr_idx(i)) = 1;
    model.A(num_constr+i, IndOfDebugVars(i)) = roundsd((-model.var_ub(Constr_idx(i))) + ConstraintTable{i,3},5);
    model.constraintType{num_constr+i} = '<';
end

% Create one more constraint involving the newly generated variable
[num_constr,~] = size(model.A);
for i=1:length(ConstraintTable(:,1))
    model.constraintNames{num_constr+i} = strcat('DebugRelax_lb_', model.varNames{Constr_idx(i)});
    model.rhs(num_constr+i) = roundsd(ConstraintTable{i,2},5);
    model.A(num_constr+i, Constr_idx(i)) = 1;
    model.A(num_constr+i, IndOfDebugVars(i)) = roundsd((-model.var_lb(Constr_idx(i))) + ConstraintTable{i,2},5);
    model.constraintType{num_constr+i} = '>';
end

indfr = find(model.f);
[num_constr, ~] = size(model.A);
model.constraintNames{num_constr+1,1} = 'Total_Violations';
model.constraintType{num_constr+1,1} = '<';
model.rhs(num_constr+1) = violations;
model.A(num_constr+1,indfr) = 1;

model.f = zeros(size(model.f,1),1);
model.f(find_cell(objrxn,model.varNames)) = 1;
model.objtype = -1;

save([NRA_dir,'model_for_alts.mat'],'model')
% Solve the system to find the minimal set of new variables that can be
% added
sol = solveTFAmodelCplex(model,timelimit);
% Get the newly added binary variables that were not able to be set to 0
if ~isempty(sol.x)
    ValuesOfDebugVars = sol.x(IndOfDebugVars);
    ids_VarToRelax = find(ValuesOfDebugVars==1);
    ids_VarToApply = find(ValuesOfDebugVars==0);
    % Those concentrations must be relaxed to get a feasible solution
    constraintsToRelax = regexprep(model.varNames(IndOfDebugVars(ids_VarToRelax)),'^DebugRelax_','');
    constraintsToApply = regexprep(model.varNames(IndOfDebugVars(ids_VarToApply)),'^DebugRelax_','');
    fc = exp(sol.val);
else
    constraintsToRelax = [];
    constraintsToApply = [];
    fc = [];
end
