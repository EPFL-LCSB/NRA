function fetchCCfromPCA(save_dir, model_name, options)

findRxn = strrep(options.objrxnPCA,'RXN_','');
if strcmp(findRxn,'Growth')==1
    load([save_dir,'MU_FCCs.mat'])
    [~, samples_out_ID,~,~,~,~] = PCA_samples(MU_FCCs);
else
    % load from previous run the CpaggDATA containing the PCA analysis
    load([save_dir,'CpaggDATA.mat'],'dataPCA_V')
    % Load kinModel from data
    load([save_dir,'metSampData',model_name,'.mat'],'kinModel')
    idRxn = find_cell(findRxn,kinModel.reactions);
    % This is the reduced coordinates in PCA space.
    Vt_samples = dataPCA_V{idRxn,1};
    
    % Center of mass coordinates which is essentially 0
    dimPCA = mean(Vt_samples,2);
    
    % Find closest point to mean of reference
    distCharVec = zeros(1,length(Vt_samples(1,:)));
    for i=1:length(Vt_samples(1,:))
        distCharVec(i) = distVec(dimPCA,Vt_samples(:,i));
    end
    [~,minValIndx] = min(distCharVec);
    refID = minValIndx;
    
    %%  Here we find extreme samples for the samples
    numPoints = 100;
    extremesID = [];
    allPC = 1:numel(dimPCA);
    for k=1:numel(dimPCA)
        %% FIND MAX PC
        [~, v_ID] = sort(Vt_samples(k,:),'descend');
        ID_samples = v_ID(1:numPoints);
        
        % Get PCA coordinates of the new representative samples so that we can
        % find the closes sample to it
        dimPC = dimPCA;
        [ID_logic_PC,~] = ismember(allPC',k);
        % Set to zero alpha so we dont move along other PCs than the one we
        % want extremes for
        dimPC = dimPC(allPC(ID_logic_PC==0));
        
        distCharVec = zeros(1,numPoints);
        for i=1:numPoints
            distCharVec(i) = distVec(dimPC,Vt_samples(allPC(ID_logic_PC==0),ID_samples(i)));
        end
        [~,minValIndx] = min(distCharVec);
        
        extremesID = [extremesID ID_samples(minValIndx)];
        
        %% FIND MIN PC
        [~, v_ID] = sort(Vt_samples(k,:),'ascend');
        ID_samples = v_ID(1:numPoints);
        
        % Get PCA coordinates of the new representative flux so that we can
        % find the closes sample to it
        dimPC = dimPCA;
        [ID_logic_PC,~] = ismember(allPC',k);
        dimPC = dimPC(allPC(ID_logic_PC==0));
        
        distCharVec = zeros(1,numPoints);
        for i=1:numPoints
            distCharVec(i) = distVec(dimPC,Vt_samples(allPC(ID_logic_PC==0),ID_samples(i)));
        end
        [~,minValIndx] = min(distCharVec);
        
        extremesID = [extremesID ID_samples(minValIndx)];
        
    end
    samples_out_ID = [refID extremesID];
end

% Now that we have the reference sample IDs and extreme IDs we can fetch
% the models from the original ORACLE files
load([save_dir,'dataInfo.mat'],'dataInfo')
load([save_dir,'metSampData',model_name,'.mat'],'kinModel')

switch options.modelcase
    case 'Ecoli'
        params = kinModel.parameters;
        reacs = params;
        params = params(cellfun(@isempty,regexp(params(:,1),'Vbs')),:);
        params = params(cellfun(@isempty,regexp(params(:,1),'LMPD')),:);
        params = params(cellfun(@isempty,regexp(params(:,1),'LCSB')),:);
        params = params(cellfun(@isempty,regexp(params(:,1),'TOe')),:);
        params = params(cellfun(@isempty,regexp(params(:,1),'_c')),:);
        params = params(cellfun(@isempty,regexp(params(:,1),'_p')),:);
        reacs = reacs(cellfun(@isempty,regexp(reacs(:,1),'Vbs')),:);
        reacs = reacs(cellfun(@isempty,regexp(reacs(:,1),'LMPD')),:);
        reacs = reacs(cellfun(@isempty,regexp(reacs(:,1),'LCSB')),:);
        reacs = reacs(cellfun(@isempty,regexp(reacs(:,1),'TOe')),:);
        reacs = reacs(cellfun(@isempty,regexp(reacs(:,1),'_c')),:);
        reacs = reacs(cellfun(@isempty,regexp(reacs(:,1),'_p')),:);
end

% location rxns and params of interest
LocationParams = find_cell(params,kinModel.parameters);
LocationRxns = find_cell(reacs,kinModel.reactions);

for i=1:numel(samples_out_ID)
    disp(['Load model no ',num2str(i), ' from ORACLE runs'])
    id_Sample = samples_out_ID(i);
    id_inBatch = mod(id_Sample,dataInfo.sizeOfBatch);
    no_Batch = floor(id_Sample/dataInfo.sizeOfBatch)+1;%if id smaller than batch size number this will be zero so we have to add 1 as it is in the first batch.
    load([save_dir,'Data',model_name,num2str(no_Batch),'.mat'],'Cpagg','Cxagg');
    % Representative CCs chosen.
    thread_FCCs = permute(Cpagg(id_inBatch,:,:),[2 3 1]);
    thread_CCCs = permute(Cxagg(id_inBatch,:,:),[2 3 1]);
    
    % reduce set CCs and theri names for preparing NRA model
    controlconc = thread_CCCs(:,LocationParams);
    controlflux = thread_FCCs(LocationRxns,LocationParams);
    consnames =  kinModel.parameters(LocationParams);
    fluxnames =  kinModel.reactions(LocationRxns);

    if options.getMu==1
        fluxnames=[fluxnames;{'Growth'}];
        load([save_dir,'MU_FCCs.mat']);
        % Fetch the growth FCC of PCA sample
        Growth_FCC=MU_FCCs(:,id_Sample);
        controlflux=[controlflux;Growth_FCC'];
    end
    
    if i==1
        aggFCC=controlflux;
        aggCCC=controlconc;
    else
        aggFCC=cat(3,aggFCC,controlflux);
        aggCCC=cat(3,aggCCC,controlconc);
    end
end

% Save the CCs for performing NRA
CCforNRA.aggFCC = aggFCC;
CCforNRA.aggCCC = aggCCC;
CCforNRA.consnames = consnames;
CCforNRA.fluxnames = fluxnames;
CCforNRA.metnames = kinModel.metabolites;
CCforNRA.samples_out_ID = samples_out_ID;

save([save_dir,'CCforNRA_',options.objrxnPCA{:},'.mat'],'CCforNRA')


%% Find distance...
    function d=distVec(a,b)
        aa=sum(a.*a,1); bb=sum(b.*b,1); ab=a'*b;
        d = sqrt(abs(aa( ones(size(bb,2),1), :)' + bb( ones(size(aa,2),1), :) - 2*ab));
    end


end