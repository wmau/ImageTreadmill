function [R,lapNum,sessionNum] = PVTrialCorr_place(mds,varargin)
%[R,lapNum,sessionNum] = PVTrialCorr_place(mds,varargin)
%
%   Performs pairwise correlations on population vectors for each trial.
%
    
%% Parse inputs.
    p = inputParser;
    p.addRequired('mds',@(x) isstruct(x)); 
    p.addParameter('useIntervals',false,@(x) islogical(x)); 
    p.addParameter('placeCellsOnly',true,@(x) islogical(x)); 
    p.addParameter('similarityMetric','corr',@(x) ischar(x));
    
    p.parse(mds,varargin{:});
    placeCellsOnly = p.Results.placeCellsOnly;
    similarityMetric = lower(p.Results.similarityMetric);
 
%% Gather data across days.
    %Get all the place cells that existed among all recording sessions.
    nSessions = length(mds); 
    cellsOfInterest = cell(nSessions,1);
    nTrials = zeros(nSessions,1);
    for s=1:nSessions
        cd(mds(s).Location);
        
        %Get place cells. 
        if placeCellsOnly
            cellsOfInterest{s} = getPlaceCells(mds(s),.01);
             %Or just all the neurons. 
        else
            load('FinalOutput.mat','NumNeurons');
            cellsOfInterest{s} = 1:NumNeurons;
        end
        
        %Get number of laps.
        load('SpatialTraces.mat','raster');
        [nTrials(s),nBins,~] = size(raster); 
        
    end
    totalTrials = sum(nTrials);
    trialBlocks = [0; cumsum(nTrials)];
    
    %Map the place cells.
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
        load('SpatialTraces.mat','raster'); 
        
        %Get neurons this session that we need to capture. 
        neurons = map(:,s);
        mapped = find(neurons>0);
        neurons = neurons(mapped);
        nMapped = length(mapped);
        
        %Capture traces.       
        for n=1:nMapped                   
            %Take the raster and transpose to make bins x trial matrix.
            cellTrace = raster(:,:,neurons(n))'; 
            PVs(:,trialBlocks(s)+1:trialBlocks(s+1),mapped(n)) = cellTrace;
        end
        nTrialsThisSession = size(raster,1);        %# trials this session.
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