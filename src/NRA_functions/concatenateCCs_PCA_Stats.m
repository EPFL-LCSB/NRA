function concatenateCCs_PCA_Stats(save_dir,my_string,model_name)
% input directory where data is located and the string for whether we
% want to concatinate cpagg (flux ccs) or cxagg (conc ccs).

dataInfo_dir=save_dir;

Cases= {model_name};

CC=my_string;% string should be Cxagg or Cpagg

% Max size of the data chunks we concatenate from the CCs.
% 20 gigs hold in my workspace kind of ok and i can still perform pca on
% this without crach. After this start becoming problematic.
maxSizeChunk=19*10^9;

for i=1:numel(Cases)
    
    load([dataInfo_dir,'metSampData',Cases{i},'.mat'],'kinModel');
    myPars=kinModel.parameters;
    
    load([dataInfo_dir,'dataInfo.mat']);
    
    if strcmp(CC,'Cpagg')
       myVars=kinModel.reactions; 
    elseif strcmp(CC,'Cxagg')
       myVars=kinModel.metabolites;
    end
        
    numVar=numel(myVars);
    numPar=numel(myPars);
    numSat=dataInfo.noRep*dataInfo.sizeOfBatch;
    
    sizeCCs=numVar*numPar*numSat*8; % estiamte how many bites the whole CC block takes
    splitRatio=sizeCCs/maxSizeChunk;
    
    numChunks=ceil(splitRatio);
    
    % Find how many variable we can take per chunk at max considering we
    % keep all parameters and samples for each chunk
    varPerChunkNum=maxSizeChunk/(numPar*numSat*8);
    
    varPerChunkInt=floor(varPerChunkNum);
    
    dataSTD=[];
    dataVAR=[];
    dataMEAN=[];
    dataMEDIAN=[];
    dataLQ=[];
    dataUQ=[];
    dataPCA=[];
    dataPCA_IDs=[];
    dataPCA_perc=[];
    dataPCA_U=[];
    dataPCA_E=[];
    dataPCA_V=[];
    for myChunk=1:numChunks
        
        if myChunk<numChunks
           el_Vars=myVars(((myChunk-1)*varPerChunkInt+1):(myChunk*varPerChunkInt));
        else
           el_Vars=myVars(((myChunk-1)*varPerChunkInt+1):end); 
        end
        
        LocationVars=find_cell(el_Vars,myVars);
        LocationPars=find_cell(myPars,myPars);
        
        % Create the Chunk matrix we populate
        allCC=zeros(dataInfo.noRep*dataInfo.sizeOfBatch,length(LocationVars),length(LocationPars));
        elastIDs=[];
        for myBatch=1:dataInfo.noRep
            load([dataInfo_dir,'Data',Cases{i},num2str(myBatch)],'elastIndAgg');
            load([dataInfo_dir,'Data',Cases{i},num2str(myBatch)],CC);
            eval(['allCC((',num2str(myBatch),'-1)*dataInfo.sizeOfBatch+1:',num2str(myBatch),...
                '*dataInfo.sizeOfBatch,:,:)=',CC,'(:,LocationVars,LocationPars);']);
            display(['Batch ',num2str(myBatch)])
            elastIDs=[elastIDs;find(elastIndAgg)+(myBatch-1)*dataInfo.sizeOfBatch];
        end
        
        allCC=allCC(elastIDs,:,:);
        
        PCAs=[];
        PCA_ids=[];
        PCAperc=[];
        PCA_U=[];
        PCA_E=[];
        PCA_V=[];
        for leVar=1:numel(el_Vars)
           CCtoPCA=allCC(:,leVar,:);
           CCtoPCA=permute(CCtoPCA,[3 1 2]);
                      
           [Samples,Samples_id,percentCover,U,E,V]=PCA_samples(CCtoPCA);

           lePCA=Samples(:,1)';
           lePCA_id=Samples_id(1);
           
           PCAs=[PCAs;lePCA];
           PCA_ids=[PCA_ids;lePCA_id];
           PCAperc=[PCAperc;percentCover];
           PCA_U=[PCA_U;{U}];
           PCA_E=[PCA_E;{E}];
           PCA_V=[PCA_V;{V}];
        end

        varRun=var(allCC,1);
        varRun=permute(varRun,[2 3 1]);
        
        stdRun=std(allCC,1);
        stdRun=permute(stdRun,[2 3 1]);
        
        meanRun=mean(allCC,1);
        meanRun=permute(meanRun,[2 3 1]);
        
        medianRun=median(allCC,1);
        medianRun=permute(medianRun,[2 3 1]);
        
        quant=quantile(allCC,0.25,1);
        lqRun=permute(quant,[2 3 1]);
        
        quant=quantile(allCC,0.75,1);
        uqRun=permute(quant,[2 3 1]);
        
        clear allCC;
        
        dataSTD=[dataSTD;stdRun];
        dataVAR=[dataVAR;varRun];
        dataMEAN=[dataMEAN;meanRun];
        dataMEDIAN=[dataMEDIAN;medianRun];
        dataLQ=[dataLQ;lqRun];
        dataUQ=[dataUQ;uqRun];
        
        dataPCA=[dataPCA;PCAs];
        dataPCA_IDs=[dataPCA_IDs;PCA_ids];
        dataPCA_perc=[dataPCA_perc;PCAperc];
        dataPCA_U=[dataPCA_U;PCA_U];
        dataPCA_E=[dataPCA_E;PCA_E];
        dataPCA_V=[dataPCA_V;PCA_V];
    end
    
    save([dataInfo_dir,CC,'DATA.mat'],'dataSTD','dataVAR','dataMEAN','dataMEDIAN','dataLQ','dataUQ',...
        'dataPCA','dataPCA_IDs','dataPCA_perc','dataPCA_U','dataPCA_E','dataPCA_V','elastIDs','-v7.3');
end

end


