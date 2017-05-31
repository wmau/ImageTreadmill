function [R,lapNum,sessionNum,Rmeans] = PVTrialCorr2(mds,varargin)
%[R,p,dayMeans] = PVTrialCorr(mds,varargin)
%
%   Performs pairwise correlations on population vectors for each treadmill
%   run. 
%
    
%% Parse inputs.
    p = inputParser;
    p.addRequired('mds',@(x) isstruct(x)); 
    p.addParameter('codingCells','timecells',@(x) ischar(x)); 
    p.addParameter('similarityMetric','corr',@(x) ischar(x));
    p.addParameter('z',true,@(x) islogical(x));
    
    p.parse(mds,varargin{:});
    codingCells = p.Results.codingCells;
    similarityMetric = lower(p.Results.similarityMetric);
    z = p.Results.z;
    
    possibleCodingCells = {'timecells','placecells','none_time','none_place'};
    assert(any(strcmp(codingCells,possibleCodingCells)),'Invalid cell type!');
 
%% Gather data across days.
    %Get all the time cells that existed among all recording sessions.
    nSessions = length(mds); 
    cellsOfInterest = cell(nSessions,1);
    nTrials = zeros(nSessions,1);
    for s=1:nSessions
        cd(mds(s).Location);
        
        switch codingCells
            case 'timecells'    %Get time cells. 
                cellsOfInterest{s} = getTimeCells(mds(s));
                
                load('TimeCells.mat','TodayTreadmillLog','T'); 
                nTrials(s) = sum(TodayTreadmillLog.complete); 
                nBins = T.*20;
            case 'placecells'   %Get place cells. 
                cellsOfInterest{s} = getPlaceCells(mds(s),0.01);
                
                load('SpatialTraces.mat','raster');
                [nTrials(s),nBins,~] = size(raster);
            case 'none_time'    %Get all cells.
                load('FinalOutput.mat','NumNeurons');
                cellsOfInterest{s} = 1:NumNeurons;    
                
                %Get number of bins and trial counts.
                load('TimeCells.mat','TodayTreadmillLog','T'); 
                nTrials(s) = sum(TodayTreadmillLog.complete); 
                
                nBins = T.*20;
            case 'none_place'   %Get all cells.
                load('FinalOutput.mat','NumNeurons');
                cellsOfInterest{s} = 1:NumNeurons;  
                
                %Get number of bins and trial counts.
                load('SpatialTraces.mat','raster');
                [nTrials(s),nBins,~] = size(raster);
        end
        
    end
    totalTrials = sum(nTrials);
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
        
        switch codingCells
            case {'timecells','none_time'}
                load('TreadmillTraces.mat','DFDTTrdmll'); 
        %         DFDTTrdmll = diff(LPtrdmll,[],2);
        %         DFDTTrdmll = padarray(DFDTTrdmll,[0 1 0],0,'pre');
                rasters = DFDTTrdmll;
            case {'placecells','none_place'}
                load('SpatialTraces.mat','raster'); 
                rasters = raster;
        end
        
        %Get neurons this session that we need to capture. 
        neurons = map(:,s);
        mapped = find(neurons>0);
        neurons = neurons(mapped);
        nMapped = length(mapped);
        
        %Capture traces. 
        for n=1:nMapped
            %Take raster and transpose to make bins x trial matrix. 
            cellTrace = rasters(:,:,neurons(n))';
            PVs(:,trialBlocks(s)+1:trialBlocks(s+1),mapped(n)) = cellTrace; 
        end
        nTrialsThisSession = size(rasters,1);
        lapNum(trialBlocks(s)+1:trialBlocks(s+1)) = 1:nTrialsThisSession;
        sessionNum(trialBlocks(s)+1:trialBlocks(s+1)) = s.*ones(1,nTrialsThisSession);
        
    end
    
    %Z-score.
    if z
        for n=1:nNeurons
            m = nanmean(PVs(:,:,n));
            sd = nanstd(PVs(:,:,n),[],1);

            PVs(:,:,n) = bsxfun(@minus,PVs(:,:,n),m);
            PVs(:,:,n) = bsxfun(@rdivide,PVs(:,:,n),sd);  
            
            temp = PVs(:,:,n);
            temp(~isfinite(temp)) = nan;
            
            PVs(:,:,n) = temp;
        end  
    end
    
%% Do correlation.
    %Preallocate.
    R = nan(totalTrials,totalTrials,nNeurons);     
    Rmeans = nan(nSessions,nSessions,nNeurons);
    switch similarityMetric
        case 'corr'
            for n=1:nNeurons
                %For each neuron, correlate trials to make TxTxN matrix.
                temp = corr(PVs(:,:,n),'rows','pairwise');
                temp(logical(eye(size(temp)))) = nan;
                
                R(:,:,n) = temp;
            end
        case 'innerproduct'
            for n=1:nNeurons
                for t1=1:totalTrials
                    for t2=t1:totalTrials
                        temp = inner(PVs(:,t1,n),PVs(:,t2,n))/nNeurons;
                        R(t1,t2,n) = temp;
                        R(t2,t1,n) = temp;
                    end
                end
                temp = R(:,:,n);
                temp(logical(eye(size(temp)))) = nan;
                
                R(:,:,n) = temp;
            end
        otherwise
            error('Invalid similarity metric!');
    end

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