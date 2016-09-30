function [R,lags] = xCorrTMRasters(md,tracetype)
%
%
%

%%
    cd(md.Location); 
    switch tracetype
        case 'FT'
            load('Pos_align.mat','FT');
        case 'rawtrace'
            load('Pos_align.mat','rawtrace');
            FT = rawtrace; clear rawtrace; 
            m = max(FT,[],2);
            FT = FT./repmat(m,1,size(FT,2)); 
        case 'difftrace'
            load('Pos_align.mat','difftrace');
            FT = difftrace; clear difftrace; 
            m = max(FT,[],2);
            FT = FT./repmat(m,1,size(FT,2));        
    end
    nNeurons = size(FT,1);
    load('TimeCells.mat','TimeCells','T','TodayTreadmillLog'); 
    
    %Make treadmill run indices even.
    inds = TodayTreadmillLog.inds; 
    inds = inds(find(TodayTreadmillLog.complete),:);
    inds(:,2) = inds(:,1) + 20*T-1; 
    
    %Get active neurons. 
    [~,neurons] = nNeuronsActiveonTM(md);
    
    %Build the rasters. 
    rasters = cell(1,nNeurons);
    for n=neurons
        switch tracetype
            case 'FT'
                rasters{n} = buildRaster(inds,FT,n);
            case {'rawtrace','difftrace'}
                rasters{n} = buildRasterTrace(inds,FT,n); 
        end
    end
    critLaps = 0.25*size(rasters{neurons(1)},1);
        
    B = 50;
    R.unsmoothed = cell(nNeurons);
    R.smoothed = cell(nNeurons);
    R.sig = cell(nNeurons);
    p = ProgressBar(length(neurons));
    if any(strcmp(tracetype,{'rawtrace','difftrace'}))
        nLaps = size(rasters{neurons(1)},1); 
    end
    
    for trig=neurons
        nonTrig = neurons(neurons > trig);
        
        for targ=nonTrig
            %Make sure the cells are active on the same lap enough. 
            if strcmp(tracetype,'FT')
                [trigLaps,~] = find(rasters{trig});
                [targLaps,~] = find(rasters{targ});
                laps = intersect(trigLaps,targLaps);
                nLaps = length(laps);  
            end

            %Do xCorr and smooth. 
            if nLaps > critLaps
                %Convert into numeric.            
                triggerRaster = single(rasters{trig}); 
                targetRaster = single(rasters{targ});
                
                switch tracetype
                    case 'FT'
                        %Trim laps. 
                        triggerRaster = triggerRaster(laps,:);
                        targetRaster = targetRaster(laps,:);
                        
                        %Cross correlation
                        [R.unsmoothed{trig,targ},lags] = xcorr_by_laps(triggerRaster,targetRaster); 
                        R.unsmoothed{trig,targ} = round(R.unsmoothed{trig,targ},3);
                        
                        R.smoothed{trig,targ} = smooth(mean(R.unsmoothed{trig,targ}))';
                    case {'rawtrace','difftrace'}
                        %Cross correlation.
                        [R.unsmoothed{trig,targ},lags] = xcorr_by_laps(rasters{trig},rasters{targ});    %Matrix for each lap.
                        R.smoothed{trig,targ} = mean(R.unsmoothed{trig,targ});                          %Trial mean. 
                        
                        %Confidence interval.                        
                        R.CI{trig,targ} = 1.96*std(R.unsmoothed{trig,targ})/sqrt(nLaps);
                end
                
                trialShuffleR = cell(1,B);
                timeShuffleR = cell(1,B); 
                trialShuffleMeans = zeros(B,length(lags));
                timeShuffleMeans = zeros(B,length(lags));
                for iter=1:B
                    %Deck-of-cards shuffle time. 
                    trigTimeShuffle = triggerRaster;
                    for l=1:nLaps
                        trigTimeShuffle(l,:) = circshift(triggerRaster(l,:),[0,randi([0,size(triggerRaster,2)])]);
                    end
                    
                    %XCorr matrix.
                    timeShuffleR{iter} = xcorr_by_laps(trigTimeShuffle,targetRaster);
                              
                    %Shuffle laps.
                    trigTrialShuff = triggerRaster(randperm(nLaps),:);
                    
                    %XCorr matrix. 
                    trialShuffleR{iter} = xcorr_by_laps(trigTrialShuff,targetRaster);
                    
                    if strcmp(tracetype,'FT')
                        trialShuffleR{iter} = round(trialShuffleR{iter},3);
                        timeShuffleR{iter} = round(timeShuffleR{iter},3);
                    end

                    %Mean curves. 
                    trialShuffleMeans(iter,:) = mean(trialShuffleR{iter});
                    timeShuffleMeans(iter,:) = mean(timeShuffleR{iter});
                end
                
                %Sort for CIs.
                trialShuffleMeans = sort(trialShuffleMeans);
                timeShuffleMeans = sort(timeShuffleMeans);
                R.trialShuffled{trig,targ}.mu = mean(trialShuffleMeans);
                R.trialShuffled{trig,targ}.upper = trialShuffleMeans(round(.95*B),:);
                R.trialShuffled{trig,targ}.lower = trialShuffleMeans(round(.05*B),:);
                
                R.timeShuffled{trig,targ}.mu = mean(timeShuffleMeans);
                R.timeShuffled{trig,targ}.upper = timeShuffleMeans(round(.95*B),:);
                R.timeShuffled{trig,targ}.lower = timeShuffleMeans(round(.05*B),:);
                
                if strcmp(tracetype,'FT')
                    R.trialShuffled{trig,targ}.mu = smooth(trialShuffleMeans{trig,targ}.mu)';
                    R.trialShuffled{trig,targ}.upper = smooth(trialShuffleMeans{trig,targ}.upper)';
                    R.trialShuffled{trig,targ}.lower = smooth(trialShuffleMeans{trig,targ}.lower)';
                end
                
                %Significant.
                R.sigTrial{trig,targ} = sum(repmat(R.smoothed{trig,targ},B,1) <= trialShuffleMeans)./B;
                R.sigTime{trig,targ} = sum(repmat(R.smoothed{trig,targ},B,1) <= timeShuffleMeans)./B;
                R.sig{trig,targ} = R.sigTrial{trig,targ} < 0.01 & R.sigTime{trig,targ} < 0.01 & ...
                    (R.smoothed{trig,targ} - R.CI{trig,targ}) > (R.trialShuffled{trig,targ}.upper) & ...
                    (R.smoothed{trig,targ} - R.CI{trig,targ}) > (R.timeShuffled{trig,targ}.upper);
                
            end     
        end
        
        p.progress;
    end
    p.stop;
    
    t = linspace(-T,T,length(lags)); 
    A = zeros(nNeurons);
    for i=neurons
        for j=neurons
            if ~isempty(R.sig{i,j})
                sigs = R.sig{i,j};
                goodt = t(sigs==1);
                goodt(goodt < -2) = nan;
                goodt(goodt > 2) = nan;
                if nanmean(goodt) < 0
                    A(i,j) = 1;
                elseif nanmean(goodt) > 0
                    A(j,i) = 1;
                end
            end
        end
    end
    
    %Divide by frame rate.
    lags = lags./20;
    
end
        