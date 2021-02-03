% Run this script in the folder of the concentration violation study you
% want to generate alternative concentration sets for the maximal objective
% value

for violations = 4 % Number of violations
    for enzymes = [2,4,7] % Number of manipulations

        load('model_for_alts.mat')
        
        maxNum = 100;
        time_limit = 3600; % time limit for the solver
        
        model.rhs(find_cell('Total_Violations',model.constraintNames)) =  violations;
        model.rhs(find_cell('TOTAL_PERTUNRATIONS',model.constraintNames)) = 271 - enzymes;
        model.constraintType{find_cell('TOTAL_PERTUNRATIONS',model.constraintNames)} = '=';

        IndOfDebugVars = getAllVar(model,{'DebugRelax'});
        save_str = ['CONC_ALTS_enz',num2str(enzymes),'.mat'];
        num_sols = 0;
        
        fprintf('Solving for alternative #1\n');
        sol = solveTFAmodelCplex(model,time_limit);
        fprintf('Fold Change = %d\n',exp(sol.val))

        model.var_lb(model.f==1) = 0.999*sol.val;
        model.var_ub(model.f==1) = 1.001*sol.val;
        
        model.f = zeros(size(model.f,1),1);
        model.f(IndOfDebugVars) = 1;

        
        while ((num_sols < maxNum) && ~(isempty(sol.x)))
            
            num_sols = num_sols + 1;
            
            ValuesOfDebugVars = sol.x(IndOfDebugVars);
            ids_VarToRelax = find(ValuesOfDebugVars==1);
            constraintsToRelax = regexprep(model.varNames(IndOfDebugVars(ids_VarToRelax)),'^DebugRelax_','');
            
            CONC_ALTS.variables(:,num_sols) = constraintsToRelax;
            cvalues = zeros(size(constraintsToRelax));
            for i=1:size(constraintsToRelax,1)
                cvalues(i) = sol.x(find_cell(constraintsToRelax(i),model.varNames));
            end
            CONC_ALTS.values(:,num_sols) = cvalues;

            [num_cons,num_vars] = size(model.A);
            new_constraint = zeros(1,num_vars);
            new_constraint(IndOfDebugVars(ids_VarToRelax)) = 1;
            model.A(num_cons+1,:) = new_constraint;
            model.rhs(num_cons+1) = violations;
            model.constraintNames{num_cons+1} = ['alt_',num2str(num_sols)];
            model.constraintType{num_cons+1} = '<';
            
            fprintf('Solving for alternative #%d\n',num_sols+1);
            sol = solveTFAmodelCplex(model,time_limit);
                        
            save(save_str,'CONC_ALTS')
        end
    end
end
