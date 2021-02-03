function getCCforMu(save_dir,model_name,options)

Cases= {model_name};

for i=1:numel(Cases)
    
    load([save_dir,'metSampData',Cases{i},'.mat'],'kinModel');
    load([save_dir,'metSampData',Cases{i},'.mat'],'tFBAmodel');
    load([save_dir,'dataInfo.mat']);
    
    %For the params we remove transfers, lumps and Vbs that aren't of interest
    reacs=kinModel.reactions;
    reacs=reacs(~cellfun(@isempty,regexp(reacs(:,1),'Vbs')),:);
    
    switch options.modelcase
        case 'Ecoli'
            params=kinModel.parameters;
            params = params(cellfun(@isempty,regexp(params(:,1),'Vbs')),:);
            params = params(cellfun(@isempty,regexp(params(:,1),'LMPD')),:);
            params = params(cellfun(@isempty,regexp(params(:,1),'LCSB')),:);
            params = params(cellfun(@isempty,regexp(params(:,1),'TOe')),:);
            params = params(cellfun(@isempty,regexp(params(:,1),'_c')),:);
            params = params(cellfun(@isempty,regexp(params(:,1),'_p')),:);
    end
    
    LocationParams = find_cell(params,kinModel.parameters);
    LocationRxns = find_cell(reacs,kinModel.reactions);
    
    allCC = zeros(dataInfo.noRep*dataInfo.sizeOfBatch,length(LocationRxns),length(LocationParams));
    
    elastIDs = [];
    for k=1:dataInfo.noRep
        load([save_dir,'Data',Cases{i},num2str(k)],'elastIndAgg');
        load([save_dir,'Data',Cases{i},num2str(k),'.mat'],'Cpagg');
        eval(['allCC((',num2str(k),'-1)*dataInfo.sizeOfBatch+1:',num2str(k),...
            '*dataInfo.sizeOfBatch,:,:)=Cpagg(:,LocationRxns,LocationParams);']);
        display(['Batch ',num2str(k)])
        elastIDs = [elastIDs;find(elastIndAgg)+(k-1)*dataInfo.sizeOfBatch];
    end
    
    allCC = allCC(elastIDs,:,:);
    
    vNet_BS = kinModel.vNet(LocationRxns);
    
    indBiomass = find_cell(options.biomassRxn,tFBAmodel.rxnsRG);
    
    vNet_mets = tFBAmodel.metShortNames(tFBAmodel.SRG(:,indBiomass)<0);
    stoich_mets = tFBAmodel.SRG(tFBAmodel.SRG(:,indBiomass)<0,indBiomass);
    
    FIND_H2O = tFBAmodel.mets(find_cell('H2O',tFBAmodel.metFormulas));
    
    StoichVnet = stoich_mets(~ismember(vNet_mets,FIND_H2O));
    
    sum_carbon = 0;
    formula = tFBAmodel.metFormulas(~ismember(vNet_mets,FIND_H2O));
    for CARBON=1:length(formula)
        try
            carbon(CARBON,1) = sscanf(formula{CARBON}, 'C%d');
        catch
            carbon(CARBON,1) = 0;
        end
        sum_carbon = sum_carbon-StoichVnet(CARBON)*carbon(CARBON,1);
    end
    Phi4BiomassByProducts = carbon/sum_carbon;
    
    Growth = tFBAmodel.SolFluxesRG(indBiomass);
    
    Construct_Mu = Phi4BiomassByProducts.*vNet_BS/Growth;
    
    MU_FCCs = zeros(length(LocationParams),numel(elastIDs));
    for k=1:numel(elastIDs)
        iCC = allCC(k,:,:);
        iCC = permute(iCC,[2 3 1]); %[rxns x params]
        MU_FCCs(:,k) = [Construct_Mu'*iCC]';% [1xrxns] * [rxns x params]
    end
    
    save([save_dir,'MU_FCCs.mat'],'MU_FCCs','-v7.3');
end





