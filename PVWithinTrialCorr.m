function [R,p] = PVWithinTrialCorr(md,varargin)
%
%
%

%%
    p = inputParser;
    p.addRequired('mds',@(x) isstruct(x));
    p.addParameter('timeCellsOnly',true,@(x) islogical(x)); 
    p.addParameter('nBins',5,@(x) isscalar(x)); 
    
    p.parse(md,varargin{:});
    timeCellsOnly = p.Results.timeCellsOnly; 
    nBins = p.Results.nBins;
    
%% 
    %Get time cells.
    cd(md.Location);
    if timeCellsOnly
        cellsOfInterest = getTimeCells(md);
    else
        load('FinalOutput.mat','NumNeurons');
        cellsOfInterest = 1:NumNeurons;
    end
    nNeurons = length(cellsOfInterest);
    
    %Get number of laps. 
    load('TimeCells.mat','TodayTreadmillLog'); 
    nLaps = sum(TodayTreadmillLog.complete); 
   
%% 
    
    load('TreadmillTraces.mat','DFDTTrdmll');
    nPointsPerTrial = size(DFDTTrdmll,2);
    binSize = nPointsPerTrial/nBins;
    chunkLims = 0:binSize:nPointsPerTrial;
    
    PVs = nan(nNeurons,nBins);
    [R,p] = deal(nan(nBins,nBins,nLaps));
    for l=1:nLaps
        
        for i=1:nBins
            interval = chunkLims(i)+1:chunkLims(i+1);
            
            for n=1:nNeurons
                cellTraces = DFDTTrdmll(l,interval,cellsOfInterest(n));

                meanActivity = mean(cellTraces); 
                PVs(n,i) = meanActivity;
            end
        end
        
        m = nanmean(PVs,2);
        sd = nanstd(PVs,[],2);
        for n=1:nNeurons
            PVs(n,:) = bsxfun(@minus,PVs(n,:),m(n));
            PVs(n,:) = bsxfun(@rdivide,PVs(n,:),sd(n));
        end
        
        [Rslice,p(:,:,l)] = corr(PVs);     
        Rslice(logical(eye(size(Rslice)))) = nan;
        R(:,:,l) = Rslice;
    end
       
end