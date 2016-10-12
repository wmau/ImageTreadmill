function [Atrl,Atrlpval,trlNullLats] = trlShfflTest(md,A,latencies)
%
%
%

%% Set up.
    cd(md.Location);
    load('TimeCells.mat','T','TodayTreadmillLog');
    load('Pos_align.mat','FT');
    
    %Trim treadmill indices.
    [inds,nRuns] = TrimTrdmllInds(TodayTreadmillLog,T);
    
    %Get roster of active neurons. 
    [~,active] = nNeuronsActiveonTM(md);
    nNeurons = size(FT,1); 
    B = 500;
    pcrit = .01;
    
    %Preallocate.
    trlNullLats = cell(nNeurons);
    Atrlpval = nan(nNeurons); 
    
%% Build rasters. 
    rasters = cell(1,nNeurons);
    for n=1:nNeurons
        rasters{n} = buildRaster(inds,FT,n);
    end
    
%% Do trial shuffle test.
    resolution = 2;
    updateInc = round(nNeurons/(100/resolution));
    p = ProgressBar(100/resolution);
    tic;
    parpool('local');
    for snk=1:nNeurons
        %If active on the treadmill...
        if ismember(snk,active)
            %Get neurons that predict it.
            el = find(A(:,snk))'; 
        
            %For each of these...
            for src=el
                %Preallocate a cell array for storing latencies derived
                %from trial randomization.
                trlShfflLats = cell(1,B); 
                
                %Find common laps.
                [snkLaps,~] = find(rasters{snk});
                [srcLaps,~] = find(rasters{src}); 
                laps = intersect(snkLaps,srcLaps); 
                nLaps = length(laps);          
                snkTemp = rasters{snk}(laps,:);
                srcTemp = rasters{src}(laps,:);
                
                %Do B trial shuffles and get latencies.
                parfor i=1:B 
                    srcShffl = srcTemp(randperm(nLaps),:);
                    trlShfflLats{i} = sjlLatFinder(srcShffl,snkTemp);
                end
                      
                %Dump into cell array. 
                trlNullLats{src,snk} = cell2mat(trlShfflLats);
                
                %Get p-value of KS test. 
                [~,Atrlpval(src,snk)] = kstest2(latencies{src,snk},...
                    trlNullLats{src,snk});
            end
        end
        
        %Update progress bar. 
        if round(snk/updateInc) == (snk/updateInc)
            p.progress;
        end

    end
    delete(gcp);
    
    Atrl = false(nNeurons); 
    Atrl(Atrlpval < pcrit) = true;
    
    elapsed = toc;
            
    save('TrialShuffle.mat','A','Atrlpval','trlNullLats','elapsed','-v7.3');
end