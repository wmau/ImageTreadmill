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
            load('TimeCells.mat','curves');
            inds = TrimTrdmllInds(TodayTreadmillLog,T);
            
            t_binned = linspace(0,10,size(curves.sig{1},2));
            t_unbinned = linspace(0,10,inds(1,2)-inds(1,1));
              
            for thisNeuron = neurons
                raster = buildRasterTrace(inds,PSAbool,thisNeuron);

                if shuffle
                    raster = raster(randperm(size(raster,1)),:);
                end
                
                sigCurve = false(1,size(raster,2));
                sigInds = find(curves.sig{thisNeuron}); 
                a = t_binned(sigInds(1));
                b = t_binned(sigInds(end)); 
                i = findclosest(t_unbinned,a); 
                j = findclosest(t_unbinned,b);
                sigCurve(i:j) = true; 
                skewness(thisNeuron) = TrialSkewness(raster,sigCurve);
            end
        case 'place'
            load('SpatialTraces','raster','sigCurve'); 
            rasters = raster; 
   
            n = 1;
            for thisNeuron = neurons
                raster = rasters(:,:,thisNeuron);

                if shuffle 
                    raster = raster(randperm(size(raster,1)),:);
                end
                skewness(thisNeuron) = TrialSkewness(raster,sigCurve(thisNeuron,:));
                n = n+1;
            end
    end
    
end