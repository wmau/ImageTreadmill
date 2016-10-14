function [Atrl,Atrlpval,trlNullLats] = trlShfflTest(md,A,latencies)
%[Atrl,Atrlpval,trlNullLats] = trlShfflTest(md,A,latencies)
%
%   Performs a trial shuffle and tests whether that produces a different
%   latency distribution between two rasters (KS test). 
%
%   INPUTS
%       md: session entry.
%
%       A: adjacency matrix. 
%
%       latencies: NxN cell array, distributions of cell to cell latencies. 
%
%   OUTPUTS
%       Atrl: adjacency matrix after trial shuffle test. 
%
%       Atrlpval: p-value of KS test of the trial shuffled latencies vs
%       real latencies.
%
%       trlNullLats: null distribution of trial shuffled latencies. 
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
    nComparisons = nNeurons * nNeurons;
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
    disp('Performing trial shuffles...');
    parpool('local');
    
    resolution = 2;
    updateInc = round(nComparisons/(100/resolution));
    p = ProgressBar(100/resolution);
    parfor c=1:nComparisons
        [src,snk] = ind2sub([nNeurons,nNeurons],c);
        
        if ismember(snk,active) && A(src,snk)
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
            for i=1:B 
                srcShffl = srcTemp(randperm(nLaps),:);
                trlShfflLats{i} = sjlLatFinder(srcShffl,snkTemp);
            end
            
            %Dump into cell array. 
            trlNullLats{c} = cell2mat(trlShfflLats);

            %Get p-value of KS test. 
            [~,Atrlpval(c)] = kstest2(latencies{c},trlNullLats{c});
        else
            Atrlpval(c) = nan;
        end        
        
        %Update progress bar. 
        if round(c/updateInc) == (c/updateInc)
            p.progress;
        end
    end
    delete(gcp);
       
    %Make adjacency matrix.
    Atrl = false(nNeurons); 
    Atrl(Atrlpval < pcrit) = true;
                
    %save('TrialShuffle.mat','Atrl','Atrlpval','trlNullLats','elapsed','-v7.3');
end