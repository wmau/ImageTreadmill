function [R,p,dayMeans] = PVTrialCorr(mds,varargin)
%
%
%

%% Parse inputs.
    p = inputParser;
    p.addRequired('mds',@(x) isstruct(x)); 
    p.addParameter('useIntervals',false,@(x) islogical(x));
    
    p.parse(mds,varargin{:});
    useIntervals = p.Results.useIntervals; 
    
%%
    nSessions = length(mds); 
    
    %Get mapping matrix. 
    mapMD = getMapMD(mds);
    load(fullfile(mapMD.Location,'batch_session_map.mat')); 
    
%% Gather data across days.
    %Get all the time cells that existed among all recording sessions.
    TCs = cell(nSessions,1);
    nTrials = zeros(nSessions,1);
    for s=1:nSessions
        %Get time cells. 
        TCs{s} = getTimeCells(mds(s));
        
        %Get number of laps.
        cd(mds(s).Location); 
        load('TimeCells.mat','TodayTreadmillLog'); 
        nTrials(s) = sum(TodayTreadmillLog.complete); 
        
    end
    totalTrials = sum(nTrials);
    
    %Map the time cells.
    map = msMatchMultiSessionCells(mds,TCs); 
    nTCs = size(map,1);
    
%% Get significance intervals. 
    sigInterval = cell(nTCs,1); 
    for s=1:nSessions      
        %Get time cells from this session.
        thisSessionTCs = getTimeCells(mds(s)); 
        
        %Get position in matrix.
        neurons = map(:,s); 
        isTC = ismember(neurons,thisSessionTCs);
        stillEmpty = cellfun('isempty',sigInterval);
        isTC = isTC & stillEmpty;
        
        cd(mds(s).Location);
        load('TreadmillTraces.mat','DFDTTrdmll');
        if useIntervals
            for n=find(isTC)'
                flat = mean(DFDTTrdmll(:,:,neurons(n)));
                sd = std(flat); 

                sigInterval{n} = find(flat>mean(flat)+2*sd);
            end
        else
            [sigInterval{:}] = deal(1:size(DFDTTrdmll,2));
        end
         
    end
    
%% Grab PVs.
    %Make big ass matrix.
    PVs = nan(nTCs,totalTrials);
    t = 1;
    for s=1:nSessions
        %Load traces. 
        cd(mds(s).Location);
        load('TreadmillTraces.mat','DFDTTrdmll'); 
        load('TimeCells.mat','curves');
        
        %Get neurons this session that we need to capture. 
        neurons = map(:,s);
        mapped = find(neurons>0);
        neurons = neurons(mapped);
        nNeurons = length(mapped);
        
        %Capture traces. 
        nLaps = size(DFDTTrdmll,1); 
        for l=1:nLaps
            for n=1:nNeurons                   
                cellTraces = DFDTTrdmll(l,sigInterval{n},neurons(n)); 
                meanActivity = mean(cellTraces);
                PVs(mapped(n),t) = meanActivity;
            
            end
            t = t+1;
        end
    end
    
    %Normalize.
    m = nanmean(PVs,2);
    sd = nanstd(PVs,[],2);
    
    for n=1:nTCs
        PVs(n,:) = bsxfun(@minus,PVs(n,:),m(n));
        PVs(n,:) = bsxfun(@rdivide,PVs(n,:),sd(n));
    end
    
%% Do correlation.
    [R,p] = corr(PVs,'rows','pairwise');
    R(logical(eye(size(R)))) = nan;
    
    figure;
    imap = imagesc(R); colormap jet; 
    set(imap,'alphadata',~isnan(R));
    hold on;
    lines = cumsum(nTrials);
    lines = [1; lines];
%     line([0 0],[0,totalTrials],'color','k','linewidth',4);
%     line([0,totalTrials],[1 1],'color','k','linewidth',4);
    for s=1:nSessions+1
        line([lines(s) lines(s)],[0,totalTrials],'color','k','linewidth',4);
        line([0,totalTrials],[lines(s) lines(s)],'color','k','linewidth',4);
    end
    axis equal; axis tight;
    
    dayMeans = zeros(nSessions,1);
    for s=1:nSessions
        thisSessionRs = R(lines(s):lines(s+1),lines(s):lines(s+1)); 
        
        dayMeans(s) = nanmean(nanmean(thisSessionRs));
    end
end