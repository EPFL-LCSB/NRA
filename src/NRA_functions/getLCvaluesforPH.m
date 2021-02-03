function [LC_lb,LC_ub,DGo_lb,DGo_ub] = getLCvaluesforPH(LC_vars,options,LC_ref,Thermo_rxns,DGo_ref)

model = options.ThermoModel;
pH = options.pHdiff;
direction = options.pHdirection;

load('DB_AlbertyUpdate_keng.mat')
thermoDB = DB_AlbertyUpdate;

model_new = model;
eval(['model_new.CompartmentData.pH = model_new.CompartmentData.pH',direction,'pH;'])
model_new = prepModelforTFA(model_new,thermoDB,model_new.CompartmentData);
model_new = convToTFA(model_new, thermoDB, [],'DGo', [], [], 0, 1, 0, 0, 0);
model_new = addNetFluxVariables(model_new);
for i=1:size(model_new.rxns,1)
    model_new.var_lb(find_cell(strcat('NF_',model_new.rxns(i)),model_new.varNames))...
        = model.var_lb(find_cell(strcat('NF_',model_new.rxns(i)),model.varNames));
    model_new.var_ub(find_cell(strcat('NF_',model_new.rxns(i)),model_new.varNames))...
        = model.var_ub(find_cell(strcat('NF_',model_new.rxns(i)),model.varNames));
end
model_new = completeMissingDGinfo(model_new, 0, 1);

LC_new = runTMinMax(model_new,LC_vars,600);

LC_lb = zeros(size(LC_vars));
LC_ub = zeros(size(LC_vars));
for i=1:size(LC_vars,1)
    LC_lb(i,1) = min(LC_new(i,1),LC_ref(i,1));
    LC_ub(i,1) = max(LC_new(i,2),LC_ref(i,2));
end


if isfield(options,'DGconstraints') && options.DGconstraints==1
    DGo_lb = zeros(size(Thermo_rxns));
    DGo_ub = zeros(size(Thermo_rxns));
    for i=1:size(Thermo_rxns,1)
        a = find_cell(strcat('DGo_', Thermo_rxns{i}),model_new.varNames);
        if isempty(a)
            a = find_cell(strcat('DGMo_', Thermo_rxns{i}),model_new.varNames);
        end
        DGo_lb(i,1) = min(model_new.var_lb(a),DGo_ref(i,1));
        DGo_ub(i,1) = max(model_new.var_ub(a),DGo_ref(i,2));
    end
end

end


