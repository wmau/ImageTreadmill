function TimefieldStats(md)
%
%
%

%%
    cd(md.Location);
    load('TimeCells.mat','curves','ratebylap','TodayTreadmillLog');
    nNeurons = length(curves.tuning);
    nLaps = sum(TodayTreadmillLog.complete);
    goodLaps = find(TodayTreadmillLog.complete)';
    
    supraThresholdRegions = cellfun(@(x) regionprops(x>0,'area','pixelidxlist'),...
        curves.tuning,'unif',0);
    [~,peakInds] = cellfun(@max,curves.tuning);
    
    regions = cell(nNeurons,1);
    for n=1:nNeurons
        nRegions = length(supraThresholdRegions{n}); 
        
        for r=1:nRegions
            if ismember(peakInds(n),supraThresholdRegions{n}(r).PixelIdxList)
                regions{n} = supraThresholdRegions{n}(r).PixelIdxList;
                continue;
            end
            
        end
    end
    
%% 
    [TFpcthits] = deal(zeros(nNeurons,1));
    for n=1:nNeurons
        [activeLap,~] = find(ratebylap(goodLaps,regions{n},n));
        TFpcthits(n) = length(unique(activeLap))/nLaps; 
    end
    
    save('TimefieldStats.mat','TFpcthits');
end