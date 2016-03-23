function [normRates,sortedRates,order] = LinearizedPFs_treadmill(MD)
%[normRate,sortedRate] = LinearizedPFs_treadmill(MD,FT)
%
%   Linearizes trajectory then computes place fields by binning FT
%   responses in space. 
%
%   INPUTS
%       MD: Session entry. 
%
%       FT: Output from Tenaspis. Can be a subset of the full neuron base. 
%
%   OUTPUTS
%       normRate: Normalized responses, non-sorted.
%
%       sortedRate: Same as normRate, but sorted by peak. 
%

%% Preliminary. 
    %Go to directory. 
    currdir = pwd; 
    cd(MD.Location); 
    
    %Get treadmill log for excluding treadmill epochs. 
    TodayTreadmillLog = getTodayTreadmillLog(MD.Animal,MD.Date,MD.Session); 
    TodayTreadmillLog = AlignTreadmilltoTracking(TodayTreadmillLog,TodayTreadmillLog.RecordStartTime);
    d = TodayTreadmillLog.direction; 
    
    %Find direction for linearizing trajectory. 
    if strfind(d,'left')
        mazetype = 'left';
    elseif strfind(d,'right')
        mazetype = 'right';
    elseif strfind(d,'alternation')
        mazetype = 'tmaze';
    end
    
    %Some parameters. 
    nBins = 60;     %Spatial bins.
    minspeed = 3;   %Speed threshold (cm/s). 
    
    %Load aligned position data. 
    load(fullfile(pwd,'Pos_align.mat'),'x_adj_cm','y_adj_cm','speed','aviFrame','FT');
    x=x_adj_cm; y=y_adj_cm; FT=logical(FT); clear x_adj_cm y_adj_cm;
    [nNeurons,nFrames] = size(FT); 
    
    %Exclude treadmill epochs. 
    inds = getTreadmillEpochs(TodayTreadmillLog,aviFrame);
    i=[];
    for e=1:size(inds,1)
        i = [i,inds(e,1):inds(e,2)];
    end
    onTM = logical(zeros(1,nFrames)); 
    onTM(i) = true; 
    
    %Speed threshold. 
    isrunning = speed>minspeed; 
    
%% Linearize trajectory and bin responses spatially.
    %Linearized trajectory. 
    X = LinearizeTrajectory_treadmill(x,y,mazetype); 
    
    %Occupancy map. 
    [occ,edges] = histcounts(X,nBins); 
    
    %Bin spatial responses. 
    rate = nan(nNeurons,nBins);
    for n=1:nNeurons
        spkpos = X(FT(n,:) & isrunning & ~onTM);
        
        binned = histcounts(spkpos,edges);
        
        rate(n,:) = binned ./ occ; 
    end
    
    %Find peak and normalize. 
    [peak,inds] = max(rate,[],2);     
    normRates = bsxfun(@rdivide,rate,peak);
      
    %Smooth. 
    sm = fspecial('gaussian'); 
    parfor n=1:nNeurons
        normRates(n,:) = imfilter(normRates(n,:),sm);
    end
    
    %Sort. 
    [~,order] = sort(inds); 
    sortedRates = normRates(order,:);
    
    cd(currdir); 
end
    