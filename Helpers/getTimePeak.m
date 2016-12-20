function [t,m] = getTimePeak(MD)
%[t,i] = getTimePeak(MD)
%
%   Loads tuning curves and searches for all the peaks. 
%
%   INPUT
%       MD entry.
%
%   OUTPUTS
%       t: Time peak. 
%
%       i: Index peak. 
%

%%
    load(fullfile(MD.Location,'TimeCells.mat'),'curves','T','ratebylap',...
        'TodayTreadmillLog');
    load(fullfile(MD.Location,'Pos_align.mat'),'FT');
    inds = TrimTrdmllInds(TodayTreadmillLog,T); 
    
    %Get tuning curve peak. 
    [~,i] = cellfun(@max,curves.tuning);
    t = i/4; 
    
    %Also get median time of response during treadmill run.
    ratebylap = ratebylap(logical(TodayTreadmillLog.complete),:,:);
    nNeurons = size(ratebylap,3);
    m = zeros(nNeurons,1);
    raster = cell(nNeurons,1);
    for n=1:nNeurons
        raster{n} = buildRaster(inds,FT,n,'onsets',false);
        [~,bin] = find(raster{n});
        
        m(n) = median(bin)/20; 
    end
    
    
end