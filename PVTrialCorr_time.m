function [R,lapNum,sessionNum] = PVTrialCorr_time(mds,varargin)
%[R,p,dayMeans] = PVTrialCorr(mds,varargin)
%
%   Performs pairwise correlations on population vectors for each treadmill
%   run. 
%
    
%% Parse inputs.
    p = inputParser;
    p.addRequired('mds',@(x) isstruct(x)); 
    p.addParameter('timeCellsOnly',true,@(x) islogical(x)); 
    p.addParameter('similarityMetric','corr',@(x) ischar(x));
    
    p.parse(mds,varargin{:});
    timeCellsOnly = p.Results.timeCellsOnly;
    similarityMetric = lower(p.Results.similarityMetric);
 
%% Gather data across days.
    %Get all the time cells that existed among all recording sessions.
    nSessions = length(mds); 
    cellsOfInterest = cell(nSessions,1);
    nTrials = zeros(nSessions,1);
    for s=1:nSessions
        cd(mds(s).Location);
        
        %Get time cells. 
        if timeCellsOnly
            cellsOfInterest{s} = getTimeCells(mds(s));
             %Or just all the neurons. 
        else
            load('FinalOutput.mat','NumNeurons');
            cellsOfInterest{s} = 1:NumNeurons;
        end
        
        %Get number of laps.
        load('TimeCells.mat','TodayTreadmillLog','T'); 
        nTrials(s) = sum(TodayTreadmillLog.complete); 
        
    end
    totalTrials = sum(nTrials);
    nBins = T.*20;
    trialBlocks = [0; cumsum(nTrials)];
    
    %Map the time cells.
    map = msMatchMultiSessionCells(mds,cellsOfInterest); 
    %map = map(all(map>0,2),:);
    nNeurons = size(map,1);
    
%% Grab PVs.
    %Make big ass matrix.
    PVs = nan(nBins,totalTrials,nNeurons);
    [lapNum,sessionNum] = deal(nan(totalTrials,1));
    for s=1:nSessions
        %Load traces. 
        cd(mds(s).Location);
        load('TreadmillTraces.mat','DFDTTrdmll'); 
%         DFDTTrdmll = diff(LPtrdmll,[],2);
%         DFDTTrdmll = padarray(DFDTTrdmll,[0 1 0],0,'pre');
        
        %Get neurons this session that we need to capture. 
        neurons = map(:,s);
        mapped = find(neurons>0);
        neurons = neurons(mapped);
        nMapped = length(mapped);
        
        %Capture traces. 
        for n=1:nMapped
            %Take raster and transpose to make bins x trial matrix. 
            cellTrace = DFDTTrdmll(:,:,neurons(n))';
            PVs(:,trialBlocks(s)+1:trialBlocks(s+1),mapped(n)) = cellTrace; 
        end
        nTrialsThisSession = size(DFDTTrdmll,1);
        lapNum(trialBlocks(s)+1:trialBlocks(s+1)) = 1:nTrialsThisSession;
        sessionNum(trialBlocks(s)+1:trialBlocks(s+1)) = s.*ones(1,nTrialsThisSession);
        
    end
    
    
%% Do correlation.
    switch similarityMetric
        case 'corr'
            %Preallocate.
            R = nan(totalTrials,totalTrials,nNeurons);      
            Rmeans = nan(nSessions,nSessions,nNeurons);
            for n=1:nNeurons
                %For each neuron, correlate trials to make TxTxN matrix.
                R(:,:,n) = corr(PVs(:,:,n),'rows','pairwise');
                
                %For each session comparison, take the mean of the
                %correlation coefficients. 
                for s1=1:nSessions
                    row = sessionNum==s1;       %Row index.
                    
                    for s2=s1:nSessions
                        col = sessionNum==s2;   %Column index.
                        
                        %Take the mean.
                        Rmeans(s1,s2,n) = nanmean(nanmean(R(row,col,n)));
                    end
                end
            end
        case 'innerproduct'
            error('Not coded yet.');
        otherwise
            error('Invalid similarity metric!');
    end
    
end