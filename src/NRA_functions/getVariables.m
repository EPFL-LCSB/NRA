function varList = getVariables(model,types)

[num_cons,num_vars] = size(model.A);

for i=1:num_vars
    temp = regexp(model.varNames(i),'_','split');
    prefix{i} = temp{1,1}{1};
end

varList = find(ismember(prefix,types))';



