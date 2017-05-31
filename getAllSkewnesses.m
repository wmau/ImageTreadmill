function skewness = getAllSkewnesses(md,varargin)
%
%
%

%% 
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x));
    p.addParameter('shuffle',false,@(x) islogical(x)); 
    
    p.parse(md,varargin{:});
    
    shuffle = p.Results.shuffle;

%%
    
    cd(md.Location);
    load('Pos_align.mat','PSAbool'); 
    load('TimeCells.mat','TodayTreadmillLog','T');
    nNeurons = size(PSAbool,1);
    
    inds = TrimTrdmllInds(TodayTreadmillLog,T);
    
    TCs = getTimeCells(md)'; 
    
    skewness = nan(nNeurons,1);
    for thisTC = TCs
        raster = buildRasterTrace(inds,PSAbool,thisTC);
        
        if shuffle
            raster = raster(randperm(size(raster,1)),:);
        end
        skewness(thisTC) = TrialSkewness(raster);
    end
    
end