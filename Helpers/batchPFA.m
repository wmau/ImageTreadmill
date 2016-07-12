function batchPFA(MDs,excluderuns)
%
%
%

%%
    nSessions = length(MDs); 
    
    for s=1:nSessions
        cd(MDs(s).Location); 
        
        excludeframes = [];
        if excluderuns          
            disp('Excluding treadmill epochs');
            try
                load(fullfile(pwd,'Pos_align.mat'),'aviFrame');
                disp('Using aviFrame from Pos_align.mat'); 
            catch
                load(fullfile(pwd,'FinalOutput.mat'),'FT');
                [~,~,~,~,~,~,aviFrame] = AlignImagingToTracking(MDs(s).Pix2CM,FT,0);
                disp('Using aviFrame from AlignImagingToTracking.');
            end
            
            load(fullfile(pwd,'TimeCells.mat'),'TodayTreadmillLog'); 
            inds = getTreadmillEpochs(TodayTreadmillLog,aviFrame);
                       
            for l=1:size(inds,1)
                excludeframes = [excludeframes, inds(l,1):inds(l,2)];
            end
        end
        
        CalculatePlacefields(MDs(s),'exclude_frames',excludeframes,'minspeed',3,...
            'alt_inputs','FinalOutput.mat');
    end
    
end