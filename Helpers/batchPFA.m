function batchPFA(MDs,excluderuns)
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
    
    for s=1:nSessions
        disp(['Analyzing ',MDs(s).Animal,' on ',MDs(s).Date,', session ',...
            num2str(MDs(s).Session),'...']);
        cd(MDs(s).Location); 
        
        %Get indices of treadmill run.
        excludeframes = [];
        if excluderuns            
            load(fullfile(pwd,'TimeCells.mat'),'TodayTreadmillLog'); 
            inds = TodayTreadmillLog.inds;
                       
            for l=1:size(inds,1)
                excludeframes = [excludeframes, inds(l,1):inds(l,2)];
            end
        end
        
        %Make place fields. 
        Placefields(MDs(s),'exclude_frames',excludeframes,...
            'Tenaspis_data','FinalOutput.mat');
    end
    
end