function [skewness,even,odd] = getAllSkewnesses(md,varargin)
%skewness = getAllSkewnesses(md,varargin)
%
%   

%% Parse inputs. 
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x));
    p.addParameter('shuffle',false,@(x) islogical(x)); 
    p.addParameter('cellType','timecells',@(x) ischar(x)); 
    p.addParameter('rasterType','time',@(x) ischar(x));
    p.addParameter('subsample',false,@(x) islogical(x)); 
    
    p.parse(md,varargin{:});
    
    shuffle = p.Results.shuffle;
    cellType = p.Results.cellType; 
    rasterType = p.Results.rasterType;
    subsample = p.Results.subsample; 

%% Load data. 
    cd(md.Location);   
    load('Pos_align.mat','PSAbool'); 
    load('TimeCells.mat','TodayTreadmillLog','T');
    nNeurons = size(PSAbool,1);
    skewness = nan(nNeurons,1);
    
%% Compute the trial skewness score for each time/place cell. 
    neurons = AcquireTimePlaceCells(md,cellType)';
    switch rasterType 
        case 'time'
            load('TimeCells.mat','curves');
            inds = TrimTrdmllInds(TodayTreadmillLog,T);         %Treadmill run indices. 
            
            t_binned = linspace(0,10,size(curves.sig{1},2));    %Vector of times, binned.
            t_unbinned = linspace(0,10,inds(1,2)-inds(1,1));    %Vector of times, unbinned. 
              
            if subsample, [even,odd] = deal(nan(nNeurons,1)); end
            %For each cell...
            for thisNeuron = neurons
                %Build a raster plot. 
                raster = buildRasterTrace(inds,PSAbool,thisNeuron);

                %Shuffle trials. 
                if shuffle, raster = raster(randperm(size(raster,1)),:); end
                
                %Complicated way of getting the statistically significant
                %tuning field. 
                sigCurve = false(1,size(raster,2));         %Preallocate.
                sigInds = find(curves.sig{thisNeuron});     %Get the indices of the curve that are significant. 
                sigStart = t_binned(sigInds(1));            %Beginning...
                sigEnd = t_binned(sigInds(end));            %...and end of the significant region.
                i = findclosest(t_unbinned,sigStart);       %Corresponding start index in unbinned vector...
                j = findclosest(t_unbinned,sigEnd);         %...and end index. 
                sigCurve(i:j) = true;                       %Mark those in the preallocated vector. 
                
                %Get skewness calculation.                      
                skewness(thisNeuron) = TrialSkewness(raster,sigCurve);
                
                %Per reviewer suggestion, take only the even/odd trials and
                %take cells that are robust, not subject to noise. 
                if subsample
                    even(thisNeuron) = TrialSkewness(raster(2:2:end,:),sigCurve); 
                    odd(thisNeuron) = TrialSkewness(raster(1:2:end,:),sigCurve); 
                end
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