function batchPFA(MDs,excluderuns,cmperbin)
%
%
%

%%
    nSessions = length(MDs); 
    
    for s=1:nSessions
        cd(MDs(s).Location); 
        
        excludeframes = [];
        if excluderuns            
            load(fullfile(pwd,'TimeCells.mat'),'TodayTreadmillLog'); 
            inds = TodayTreadmillLog.inds;
                       
            for l=1:size(inds,1)
                excludeframes = [excludeframes, inds(l,1):inds(l,2)];
            end
        end
        
        CalculatePlacefields(MDs(s),'exclude_frames',excludeframes,'minspeed',3,...
            'alt_inputs','FinalOutput.mat','cmperbin',cmperbin,'HalfWindow',0);
    end
    
end