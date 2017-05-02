function [binnedCoeffMeans,binnedCoeffs] = binCoeffs(R,varargin)
%[binnedCoeffMeans,binnedCoeffs] = binCoeffs(R,varargin)
%
%   Bin correlation coefficients in R in a way specified by processingMode.
%   
%   binByDay...requires R and sessionNum.
%       R is a TxT matrix where T is the number of trials across all the
%       days.
%       sessionNum is a vector containing session number identity. 
%
%   binByNTrials...requires R and lapNum.
%       R is a txt matrix where t is the number of trials in ONE day. 
%       lapNum is a vector containing lap number identity. 
%
%   OUTPUT
%       binnedCoeffs: BxB (B=number of specified bins or days, usually 5)
%       matrix containing the mean correlation coefficient for trial blocks
%       or days. 
%

%% Input processing.
    p = inputParser;
    p.addRequired('R',@(x) isnumeric(x)); 
    p.addParameter('lapNum',[],@(x) isnumeric(x)); 
    p.addParameter('sessionNum',[],@(x) isnumeric(x)); 
    p.addParameter('nBins',5,@(x) isscalar(x)); 
    p.addParameter('processingMode','none',@(x) ischar(x)); 
    
    p.parse(R,varargin{:});
    
    lapNum = p.Results.lapNum;
    sessionNum = p.Results.sessionNum;
    nBins = p.Results.nBins;
    processingMode = p.Results.processingMode;
    
    %processingMode must be either these two:
    possibleModes = {'binByDay','binByNTrials'};
    assert(any(strcmp(processingMode,possibleModes)),...
        'Error: processingMode must be either ''binByDay'' or ''binByNTrials''.');
    
%% 
    binnedCoeffMeans = nan(nBins);
    binnedCoeffs = cell(nBins);
    switch processingMode
        case 'binByDay'
            %Error check.
            assert(~isempty(sessionNum),'Error: Needed session numbers.');
            
            %Get number of sessions.
            nSessions = max(unique(sessionNum));
            
            for s1=1:nSessions
                for s2=s1:nSessions              
                    %Grab the chunk.
                    row = sessionNum==s1;
                    col = sessionNum==s2;
                    chunk = R(row,col); 
                    
                    %Only get upper triangle. 
                    inds = logical(triu(ones(size(chunk))));
                    chunk = chunk(inds);
                    flatChunk = chunk(:);
                    
                    %Append onto diagonal.
                    binnedCoeffs{s1,s2} = [binnedCoeffs{s1,s2}; flatChunk];          
                    if s2~=s1
                        binnedCoeffs{s2,s1} = [binnedCoeffs{s2,s1}; flatChunk];
                    end
                    
                    %Take the mean.
                    m = nanmean(flatChunk);
                    binnedCoeffMeans(s1,s2) = m;
                    binnedCoeffMeans(s2,s1) = m;
                end
            end
            
        case 'binByNTrials'
            %Error check.
            assert(~isempty(lapNum),'Error: Needed lap numbers.');
            
            %Get number of laps. 
            nLaps = max(lapNum);
            
            %Make vector describing bin limits.
            binSize = round(nLaps/(nBins+1));
            chunkLims = 0:binSize:nLaps;
            chunkLims(end) = nLaps;
            
            for b1=1:nBins
                for b2=b1:nBins
                    %Grab the chunk.
                    row = chunkLims(b1)+1:chunkLims(b1+1);
                    col = chunkLims(b2)+1:chunkLims(b2+1); 
                    chunk = R(row,col);
                    
                    %Only get upper triangle. 
                    inds = logical(triu(ones(size(chunk))));
                    chunk = chunk(inds);
                    flatChunk = chunk(:);
                    
                    %Take the mean.
                    m = nanmean(flatChunk);
                    binnedCoeffMeans(b1,b2) = m;
                    binnedCoeffMeans(b2,b1) = m;
                end
            end
    end
                
end