function tempInfo(MD)
%[I,sig] = tempInfo(MD)
%
%   Calculates the Shannon mutual information I(X,K) between the random
%   variables spike count [0,1] and time via the equations: 
%
%   (1) I_time(ti) = sum[k>=0](P_k|ti * log2(P_k|ti / P_k)) 
%
%   (2) MI = sum[i=1->N](P_xi * I_time(ti)
%
%   where
%       P_ti is the probability the mouse is in time bin ti,
%       (1/(T*20)).*ones(1,T*20)
%       
%       P_k is the probability of observing k spikes,
%       sum(rasters{n}(:))./ prod(size(rasters{n})
%
%       P_k|ti is the conditional probability of observing k spikes at time
%       ti, mean(rasters{n})
%


%% Load.
    cd(MD.Location); 
    load(fullfile(pwd,'TimeCells.mat'),'TodayTreadmillLog','T');
    load(fullfile(pwd,'Pos_align.mat'),'PSAbool');
    nNeurons = size(PSAbool,1);
    B = 1000;
    totalLaps = sum(TodayTreadmillLog.complete);
    prop = 0.25;
    crit = round(prop*totalLaps);
    
%% Build spike raster.    
    %Build binned spike raster. 
    inds = TrimTrdmllInds(TodayTreadmillLog,T);
    rasters = cell(nNeurons,1);
    for n=1:nNeurons
        rasters{n} = buildRaster(inds,PSAbool,n,'onsets',false);
    end
    
%% Compute temporal information and permutation test. 
    MI = zeros(nNeurons,1);             %Preallocate empirical temporal information [bits/sec].
    Isec = zeros(nNeurons,1);           %Preallocate firing rate. 
    Ispk = zeros(nNeurons,1);           %Preallocate specificity [bits/spike].
    
    resolution = 2;
    updateInc = round(nNeurons/(100/resolution));
    p = ProgressBar(100/resolution);
    parfor n=1:nNeurons
        %Empirical temporal information. 
        [MI(n),Isec(n),Ispk(n)] = tempInfoOneNeuron(rasters{n});
        
        for i=1:B
            %Card shuffle. 
            ratebylapShift = permuteTime(rasters{n});
            
            %Compute temporal information of permuted ratebylap. 
            surrogate(n,i) = tempInfoOneNeuron(ratebylapShift);
        end
    
        %Update progress bar. 
        if round(n/updateInc) == (n/updateInc)
            p.progress;
        end
    end
    p.stop;
    
    %P-value of temporal information being significantly better than
    %shuffled data. 
    pval = sum(surrogate>=repmat(MI,[1,B]),2)/B;
    [nLapsActive,~] = cellfun(@find,rasters,'unif',0);
    nLapsActive = cellfun('length',cellfun(@unique,nLapsActive,'unif',0));
    
    %Significant neurons. 
    sig = pval<0.01 & MI>0 & nLapsActive>crit;
    
    save('TemporalInfo.mat','MI','Isec','Ispk','pval','sig');
end