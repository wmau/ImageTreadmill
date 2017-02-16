function [rate,normRates,sortedRates,order,X,edges,peakInds] = LinearizedPFs_treadmill(MD,varargin)
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

%% Parse inputs.
    p = inputParser;
    p.addRequired('MD',@(x) isstruct(x));
    p.addParameter('plotit',true,@(x) islogical(x));
    p.addParameter('PlaceCells',getPlaceCells(MD,.01),@(x) isnumeric(x));
    p.addParameter('order',false);
    
    p.parse(MD,varargin{:});
    
    plotit = p.Results.plotit;
    PCs = p.Results.PlaceCells;
    order = p.Results.order; 
    
%% Preliminary. 
    %Go to directory. 
    currdir = pwd; 
    cd(MD.Location); 
    
    %Get treadmill log for excluding treadmill epochs. 
    load('TimeCells.mat','TodayTreadmillLog'); 
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
    nBins = 80;     %Spatial bins.
    minspeed = 3;   %Speed threshold (cm/s). 
    
    %Load aligned position data. 
    load(fullfile(pwd,'Pos_align.mat'),'x_adj_cm','y_adj_cm','speed','time_interp','PSAbool');
    x=x_adj_cm; y=y_adj_cm; PSAbool=logical(PSAbool); clear x_adj_cm y_adj_cm;
    [nNeurons,nFrames] = size(PSAbool); 
    
    %Exclude treadmill epochs. 
    peakInds = TodayTreadmillLog.inds;
    i=[];
    for e=1:size(peakInds,1)
        i = [i,peakInds(e,1):peakInds(e,2)];
    end
    onTM = false(1,nFrames); 
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
        spkpos = X(PSAbool(n,:) & isrunning & ~onTM);
        
        binned = histcounts(spkpos,edges);
        
        rate(n,:) = binned ./ occ; 
    end
      
    %Find peak and normalize. 
    [peak,peakInds] = max(rate,[],2);     
    normRates = bsxfun(@rdivide,rate,peak);
      
    %Smooth. 
    sm = fspecial('gaussian'); 
    parfor n=1:nNeurons
        normRates(n,:) = imfilter(normRates(n,:),sm);
    end
    
    %Sort. 
    if ~order, [peakInds,order] = sort(peakInds(PCs)); 
    else, temp = peakInds(PCs); peakInds = temp(order); end
    sortedRates = normRates(PCs(order),:);
    
    if plotit
        imagesc(1:size(sortedRates,2),1:length(PCs),sortedRates);
        set(gca,'ydir','reverse'); axis tight;
        colormap hot;
        xlabel('Linearized Distance');
        ylabel('Neurons');
    end
    
    cd(currdir); 
end
    