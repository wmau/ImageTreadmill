function batchPFA_onTM(MDs)
%batchPFA(MDs,excluderuns,cmperbin)
%
%   Performs batch place field analysis on multiple sessions. Wrapper
%   function for Placefields. Excludes epochs of treadmill running. 
%
%   INPUTS
%       MDs: sessions to analyze.
%
%       excluderuns: logical for excluding treadmill run epochs.
%
%       cmperbin: I normally use 1. 
%

%% For each session, get indices of treadmill running and run Placefields. 
    nSessions = length(MDs); 
    nSeconds = 0;
    extraExclude = 20*nSeconds;
    
    for s=1:nSessions
        disp(['Analyzing ',MDs(s).Animal,' on ',MDs(s).Date,', session ',...
            num2str(MDs(s).Session),'...']);
        cd(MDs(s).Location); 
        
        %Get indices of treadmill run.
        excludeframes = [];
        load(fullfile(pwd,'TimeCells.mat'),'TodayTreadmillLog'); 
        inds = TodayTreadmillLog.inds;

        load('Pos_align.mat','PSAbool');
        nFrames = size(PSAbool,2);

        for l=1:size(inds,1)
            window = (inds(l,1)-extraExclude):(inds(l,2)+extraExclude);
            window = window(window>1);
            window = window(window<nFrames);

            excludeframes = [excludeframes, window];
        end
        excludeframes = setdiff(1:nFrames,excludeframes); 
        
        %Make place fields. 
        Placefields(MDs(s),'exclude_frames',excludeframes,...
            'Tenaspis_data','FinalOutput.mat','cmperbin',1,'saveName',...
            'Placefields_onTM.mat','saveSI',false,'B',1,'minspeed',0);
    end
    
end