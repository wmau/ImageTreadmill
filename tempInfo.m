function [I,sig] = tempInfo(MD)
%[I,sig] = tempInfo(MD)
%
%   Calculates the temporal information of all the neurons in a session. 
%
%   INPUT
%       MD: Session entry. 
%
%   OUPUTS
%       I: Temporal information.
%
%       sig: Significantly high temporal information, compared to permuted
%       data. 
%

%% Set up.
    currdir = pwd; 
    cd(MD.Location); 
    load(fullfile(pwd,'TimeCells.mat'),'TodayTreadmillLog','T');
    load(fullfile(pwd,'Pos_align.mat'),'FT');
    nNeurons = size(FT,1);
    
    %Number of shuffles. 
    B = 1000;
    
    %Build binned spike raster. 
    inds = TrimTrdmllInds(TodayTreadmillLog,T);
    rasters = cell(nNeurons,1);
    ratebylap = cell(nNeurons,1);
    edges = [0:20:T*20];
    for n=1:nNeurons
        rasters{n} = buildRaster(inds,FT,n,'sprs',false);
        ratebylap{n} = zeros(size(rasters{n},1),length(edges)-1);
        for l=1:size(rasters{n},1)
            
            ts = find(rasters{n}(l,:));         
            ratebylap{n}(l,:) = histcounts(ts,edges);
        end
    end
    
%% Compute temporal information and permutation test. 
    I = zeros(nNeurons,1);              %Preallocate empirical temporal information.
    surrogate = zeros(nNeurons,B);      %Preallocate surrogate data. 
    
    %parpool('local');
    resolution = 2;
    updateInc = round(nNeurons/(100/resolution));
    p = ProgressBar(100/resolution);
    parfor n=1:nNeurons
        %Empirical temporal information. 
        I(n) = tempInfoOneNeuron(ratebylap{n});
        
        for i=1:B
            %Card shuffle. 
            ratebylapShift = permuteTime(ratebylap{n});
            
            %Compute temporal information of permuted ratebylap. 
            surrogate(n,i) = tempInfoOneNeuron(ratebylapShift);
        end
    
        %Update progress bar. 
        if round(n/updateInc) == (n/updateInc)
            p.progress;
        end
    end
    p.stop;
    %delete(gcp);
    
    %P-value of temporal information being significantly better than
    %shuffled data. 
    pval = sum(surrogate>repmat(I,[1,B]),2)/B;
    
    %Significant neurons. 
    sig = pval<0.05 & I>0;
    
    save('TemporalInfo.mat','I','sig');
    cd(currdir); 
end

function I = tempInfoOneNeuron(raster)
    [nLaps,nBins,~] = size(raster);
    
    p = 1/nBins;                                    %Probability occupancy in time (uniform).
    lambda_i = mean(raster);                   %Mean rate at this time bin. (multiply by 4 because 0.25 time bins)
    lambda = sum(sum(raster))/(nLaps*nBins);   %Mean rate. 
    
    %Temporal information.
    I = nansum(p.*lambda_i.*log2(lambda_i./lambda));
end