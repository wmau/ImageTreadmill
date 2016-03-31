function treadmillRhythmicity(MD)
%
%
%

%% 
    cd(MD.Location); 
    load(fullfile(pwd,'Pos_align.mat'),'FT','aviFrame'); 
    nFrames = length(aviFrame);
    
    TodayTreadmillLog = getTodayTreadmillLog(MD.Animal,MD.Date,MD.Session); 
    TodayTreadmillLog = AlignTreadmilltoTracking(TodayTreadmillLog,TodayTreadmillLog.RecordStartTime);

    running = getTreadmillEpochs(TodayTreadmillLog,aviFrame); 
    nLaps = size(running,1); 
    L = mode(diff(running,[],2)); 
    ts = cell(nLaps,1);
    for l=1:nLaps
        ts{l} = running(l,1):running(l,2); 
        
        %Even out the lengths. 
        if length(ts{l}) > L
            d = length(ts{l}) - L - 1; 
            ts{l}(end-d:end) = [];
        elseif length(ts{l}) < L
            ts{l}(end+1) = ts{l}(end)+1;
        end
    end
    
    %Get completed laps only. 
    goodLaps = find(TodayTreadmillLog.complete)';
    nGood = length(goodLaps);
    
    P1 = zeros(nGood,L/2+1);
    for l=goodLaps
        [~,P1(l,:)] = synchFFT(FT,ts{l});
    end
    
    onTM = cell2mat(ts(goodLaps));
    sample = 1:nFrames-L; 
    sample(ismember(sample,onTM(:))) = []; 
    B = 1000; 
    rs = randsample(sample,B);
    for r=1:B
        [f,Pr(r,:)] = synchFFT(FT,rs(r):rs(r)+L-1);
    end
    
    keyboard;
end