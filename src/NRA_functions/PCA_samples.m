function [samples_out, samples_out_ID,percentCover,U,E,V] = PCA_samples(my_samples,PC_cutoff,numPoints)
% [samples_out, samples_out_ID] = PCA_samples(my_samples,PC_cutoff,numPoints)
% my_samples - [m x n] matrix with variables x number of samples 
% PC_cutoff  - Number of prinipal components used to pick extremes.
%              Function wil output two extremes per principal component and the
%              reference case. Default is 9.
% numPoints  - number of top samples at each extreme considered for extreme
%              Default is 100.

% samples_out - [reference_sample max_PC1_sample min_PC1_sample max_PC2_sample min_PC2_sample ....]
% samples_out_ID - [reference_ID max_PC1_ID min_PC1_ID max_PC2_ID min_PC2_ID ....]

if nargin<3
    numPoints = 100;
end

if nargin<2
    PC_cutoff = 9;
end

if numPoints>size(my_samples,2)
   error('You trying to use more points to pick extremes than you have samples!'); 
end

%% Reference case for samples

% Subtract mean of samples from each sample to centre the data
dataSetMean = mean(my_samples,2);
dataSetCenter = my_samples-repmat(dataSetMean,1,(size(my_samples,2)));

% Perform singular value decomposition
[U,E,V] = svd(dataSetCenter,'econ');
sigmas = diag(E);
varContrib = 100*sigmas./sum(sigmas);
cutOffTh = PC_cutoff;
percentCover = sum(varContrib(1:cutOffTh));
fprintf ('Contribution of the first %d components (percentwise): %g \n',cutOffTh,percentCover)
fprintf('Contribution of the %d th component: %g \n',cutOffTh+1,varContrib(cutOffTh+1))

% Extract the principal components
modes = U(:,1:cutOffTh);
% Compute the Loading1s
Y = modes'*dataSetCenter;
% Generate estimates from a reduced set of components
Xtest = modes*Y;

% Transpose V
Vt_samples = V';

% Weigths for mean point
dimPCA = mean(Vt_samples(1:PC_cutoff,:),2);

% Find closest point to mean of reference
distCharVec = zeros(1,length(Vt_samples(1,:)));
for i=1:length(Vt_samples(1,:))
    distCharVec(i) = distVec(dimPCA,Vt_samples([1:PC_cutoff],i));
end
[~,minValIndx] = min(distCharVec);
refSamp = my_samples(:,minValIndx);
refID = minValIndx;

%%  Here we find extreme samples for the samples
extremesSamp = [];
extremesID = [];
allPC = 1:PC_cutoff;
for k=1:PC_cutoff  
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
    dataSet_match = my_samples(:,ID_samples(minValIndx));
    extremesSamp = [extremesSamp dataSet_match];
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
    dataSet_match = my_samples(:,ID_samples(minValIndx));
    extremesSamp = [extremesSamp dataSet_match];
    extremesID = [extremesID ID_samples(minValIndx)];
    
end

samples_out = [refSamp extremesSamp];
samples_out_ID = [refID extremesID];

end

%% Find distance...
function d=distVec(a,b)
aa = sum(a.*a,1); bb = sum(b.*b,1); ab = a'*b;
d = sqrt(abs(aa( ones(size(bb,2),1), :)' + bb( ones(size(aa,2),1), :) - 2*ab));
end