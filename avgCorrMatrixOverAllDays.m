function binnedCoeffs = avgCorrMatrixOverAllDays(R,lapNum,sessionNum,processingMode,varargin)
%binnedCoeffs = avgCorrMatrixOverAllDays(R,lapNum,sessionNum,processingMode)
%
%   Takes the average across multiple binned correlation matrices. The next
%   step in the pipeline after PVTrialCorr.
%
%   INPUTS
%       R: Giant correlation matrix for one animal containing all the
%       comparisons between all sessions. 
%
%       lapNum: Vector the same length as R containing lap number
%       identities. 
%
%       sessionNum: Same as lapNum but for session numbers.
%
%       processingMode: Either 'binByNTrials' or 'binByDay'. 
%           binByNTrials: Bins the trials on one day into 5 equal chunks
%           then averages them via binCoeffs. This function then takes the
%           grand average over all sessions. 
%   
%           binByDay: Bins the trials into 5 days then averages them via
%           binCoeffs. Take the grand average after output.
%

%% Setup.
    p = inputParser;
    p.addRequired('R',@(x) isnumeric(x)); 
    p.addRequired('lapNum',@(x) isnumeric(x)); 
    p.addRequired('sessionNum',@(x) isnumeric(x)); 
    p.addRequired('processingMode',@(x) ischar(x));
    p.addParameter('nBins',5,@(x) isscalar(x)); 
    
    p.parse(R,lapNum,sessionNum,processingMode,varargin{:});
    
    nBins = p.Results.nBins;
    
    nSessions = max(sessionNum);
    binnedCoeffs = nan(nBins,nBins,nSessions);
    
%% Average over sessions.
    switch processingMode
        case 'binByNTrials'
            for s=1:nSessions
                %Get only the laps in the current session.
                theseLaps = sessionNum==s;
                laps = lapNum(theseLaps);

                %Get only that chunk of R. 
                chunk = R(theseLaps,theseLaps);

                %Bin and average.
                binnedCoeffs(:,:,s) = binCoeffs(chunk,'processingMode','binByNTrials',...
                    'lapNum',laps,'nBins',nBins);
            end    
        case 'binByDay'
            %Bin by day.
            binnedCoeffs = binCoeffs(R,'processingMode','binByDay','sessionNum',...
                sessionNum);
        otherwise 
            error('Wrong processing mode!');
    end
end