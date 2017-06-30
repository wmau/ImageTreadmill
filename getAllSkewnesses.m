function skewness = getAllSkewnesses(md,varargin)
%
%
%

%% 
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x));
    p.addParameter('shuffle',false,@(x) islogical(x)); 
    p.addParameter('cellType','timecells',@(x) ischar(x)); 
    p.addParameter('rasterType','time',@(x) ischar(x)); 
    
    p.parse(md,varargin{:});
    
    shuffle = p.Results.shuffle;
    cellType = p.Results.cellType; 
    rasterType = p.Results.rasterType;

%%
    cd(md.Location);   
    load('Pos_align.mat','PSAbool'); 
    load('TimeCells.mat','TodayTreadmillLog','T');
    nNeurons = size(PSAbool,1);
    skewness = nan(nNeurons,1);
    
    neurons = AcquireTimePlaceCells(md,cellType)';
    switch rasterType 
        case 'time'
            inds = TrimTrdmllInds(TodayTreadmillLog,T);
              
            for thisNeuron = neurons
                raster = buildRasterTrace(inds,PSAbool,thisNeuron);

                if shuffle
                    raster = raster(randperm(size(raster,1)),:);
                end
                skewness(thisNeuron) = TrialSkewness(raster);
            end
        case 'place'
            load('SpatialTraces','raster'); 
            rasters = raster; 
   
            n = 1;
            for thisNeuron = neurons
                raster = rasters(:,:,thisNeuron);

                if shuffle 
                    raster = raster(randperm(size(raster,1)),:);
                end
                skewness(thisNeuron) = TrialSkewness(raster);
                n = n+1;
            end
    end
    
end