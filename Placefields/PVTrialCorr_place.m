function [R,lapNum,sessionNum] = PVTrialCorr_place(mds,varargin)
%[R,lapNum,sessionNum] = PVTrialCorr_place(mds,varargin)
%
%   Performs pairwise correlations on population vectors for each treadmill
%   run. 
%
    
%% Parse inputs.
    p = inputParser;
    p.addRequired('mds',@(x) isstruct(x)); 
    p.addParameter('useIntervals',false,@(x) islogical(x)); 
    p.addParameter('placeCellsOnly',true,@(x) islogical(x)); 
    p.addParameter('plotit',false,@(x) islogical(x));
    p.addParameter('similarityMetric','corr',@(x) ischar(x));
    
    p.parse(mds,varargin{:});
    placeCellsOnly = p.Results.placeCellsOnly;
    plotit = p.Results.plotit;
    similarityMetric = lower(p.Results.similarityMetric);
 
%% Gather data across days.
    %Get all the time cells that existed among all recording sessions.
    nSessions = length(mds); 
    cellsOfInterest = cell(nSessions,1);
    nTrials = zeros(nSessions,1);
    for s=1:nSessions
        cd(mds(s).Location);
        
        %Get time cells. 
        if placeCellsOnly
            cellsOfInterest{s} = getPlaceCells(mds(s),.01);
             %Or just all the neurons. 
        else
            load('FinalOutput.mat','NumNeurons');
            cellsOfInterest{s} = 1:NumNeurons;
        end
        
        %Get number of laps.
        load('SpatialTraces.mat','raster');
        nTrials(s) = size(raster,1); 
        
    end
    totalTrials = sum(nTrials);
    
    %Map the time cells.
    map = msMatchMultiSessionCells(mds,cellsOfInterest); 
    %map = map(all(map>0,2),:);
    nNeurons = size(map,1);
    
%% Grab PVs.
    %Make big ass matrix.
    PVs = nan(nNeurons,totalTrials);
    [lapNum,sessionNum] = deal(nan(totalTrials,1));
    t = 1;
    for s=1:nSessions
        %Load traces. 
        cd(mds(s).Location);
        load('SpatialTraces.mat','raster'); 
        
        %Get neurons this session that we need to capture. 
        neurons = map(:,s);
        mapped = find(neurons>0);
        neurons = neurons(mapped);
        nMapped = length(mapped);
        
        %Capture traces. 
        nLaps = size(raster,1); 
        for l=1:nLaps
            for n=1:nMapped                   
                cellTrace = raster(l,:,neurons(n)); 
                meanActivity = nanmean(cellTrace);
                PVs(mapped(n),t) = meanActivity;
            
            end
            lapNum(t) = l;
            sessionNum(t) = s;
            
            t = t+1;  
        end
    end
    
    %Normalize.
    m = nanmean(PVs,2);
    sd = nanstd(PVs,[],2);
    
    for n=1:nNeurons
        PVs(n,:) = bsxfun(@minus,PVs(n,:),m(n));
        PVs(n,:) = bsxfun(@rdivide,PVs(n,:),sd(n));
    end
    
%% Do correlation.
    switch similarityMetric
        case 'corr'
            R = corr(PVs,'rows','pairwise');
        case 'innerproduct'
            R = nan(totalTrials);
            for t1=1:totalTrials
                for t2=1:totalTrials
                    R(t1,t2) = inner(PVs(:,t1),PVs(:,t2))/nNeurons;
                end
            end
        otherwise
            error('Invalid similarity metric!');
    end
    R(logical(eye(size(R)))) = nan;
    
    if plotit
        figure;
        imap = imagesc(R); colormap jet; 
        set(imap,'alphadata',~isnan(R));
        hold on;
        newDayStartInd = cumsum(nTrials);
        newDayStartInd = [1; newDayStartInd];
    %     line([0 0],[0,totalTrials],'color','k','linewidth',4);
    %     line([0,totalTrials],[1 1],'color','k','linewidth',4);
        for s=1:nSessions+1
            line([newDayStartInd(s) newDayStartInd(s)],[0,totalTrials],'color','k','linewidth',4);
            line([0,totalTrials],[newDayStartInd(s) newDayStartInd(s)],'color','k','linewidth',4);
        end
        axis equal; axis tight;
    end
    
end