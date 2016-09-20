function [X,y] = SourceSinkGLMSetUp(md,sources,sink,tracetype)
%[X,y] = SourceSinkGLMSetUp(md,sources,sink)
%
%   Sets up design matrix and response variable given neuron indices. 
%
%   INPUTS
%       md: Session entry.
%
%       sources (optional): Vector of neuron indices whose vectorized
%       rasters will go into X. 
%
%       sink: Neuron index whose vectorized raster will go into Y. 
%
%   OUTPUTS
%       X: Design matrix of predictor variables. First column is time,
%       columns after that are vectorized rasters. 
%
%       y: Response variable, binary. Vectorized raster of sink. 
%

%% Set up.  
        cd(md.Location);
        if strcmp(tracetype,'FT')   
            load('Pos_align.mat','FT'); 
        elseif strcmp(tracetype,'trace')
            load('Pos_align.mat','trace');
        end
        load('TimeCells.mat','TimeCells','T','TodayTreadmillLog'); 

        %Make treadmill run indices even.
        inds = TodayTreadmillLog.inds; 
        inds = inds(find(TodayTreadmillLog.complete),:);
        inds(:,2) = inds(:,1) + 20*T-1; 

        if strcmp(tracetype,'FT')
            sinkRaster = buildRaster(inds,FT,sink);
        elseif strcmp(tracetype,'trace')
            sinkRaster = buildRasterTrace(inds,trace,sink);
        end
        
    if nargout > 1
        [nLaps,nTimeBins] = size(sinkRaster); 
        t = linspace(0,T,nTimeBins)';
        X(:,1) = repmat(t,nLaps,1);

        nSources = length(sources);
        for s=1:nSources
            if strcmp(tracetype,'FT')
                sourceRaster = buildRaster(inds,FT,sources(s));
            elseif strcmp(tracetype,'trace')
                sourceRaster = buildRasterTrace(inds,trace,sources(s));
            end
            sourceRaster = sourceRaster';
            X(:,s+1) = sourceRaster(:);
        end
        
        tbltitles = cell(1,nSources+1);
        tbltitles{1} = 't';
        
        for s=1:nSources
            tbltitles{s+1} = ['n',num2str(sources(s))];
        end
        
        X = array2table(X,'variablenames',tbltitles);
        X.Properties.VariableDescriptions{1} = 't';
        
        for s=1:nSources
            X.Properties.VariableDescriptions{s+1} = num2str(sources(s));
        end
    end
        
    y = sinkRaster';
    y = y(:);
    
end