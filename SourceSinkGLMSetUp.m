function [X,y] = SourceSinkGLMSetUp(md,sources,sink,lags,tracetype)
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
            load('Pos_align.mat','FT');         %Load FT.
        elseif strcmp(tracetype,'rawtrace')
            load('Pos_align.mat','rawtrace');   %Load rawtrace.
            
            %Normalize to max. 
            m = max(rawtrace,[],2);
            rawtrace = rawtrace./repmat(m,1,size(rawtrace,2));
        end
        load('TimeCells.mat','TimeCells','T','TodayTreadmillLog'); 

        %Make treadmill run indices even.
        inds = TodayTreadmillLog.inds; 
        inds = inds(find(TodayTreadmillLog.complete),:);
        inds(:,2) = inds(:,1) + 20*T-1;

        %Build raster depending on whether the input is FT or rawtrace. 
        if strcmp(tracetype,'FT')
            sinkRaster = buildRaster(inds,FT,sink);
        elseif strcmp(tracetype,'rawtrace')
            sinkRaster = buildRasterTrace(inds,rawtrace,sink);
        end
        
    %If more than one output is specified, assume user wants X and y. 
    if nargout > 1
        [nLaps,nTimeBins] = size(sinkRaster); 
        
        %Make time vector for each lap. 
        t = linspace(0,T,nTimeBins)';
        X(:,1) = repmat(t,nLaps,1);

        %For each source. 
        nSources = length(sources);
        nLags = length(lags);
        n = 2;
        for s=1:nSources
            for l=1:nLags
                if strcmp(tracetype,'FT')
                    sourceRaster = buildRaster(inds+lags(l),FT,sources(s));
                elseif strcmp(tracetype,'rawtrace')
                    sourceRaster = buildRasterTrace(inds+lags(l),rawtrace,sources(s));
                end
                sourceRaster = sourceRaster';
                X(:,n) = sourceRaster(:);
                
                n = n+1;
            end
        end
        
        tbltitles = cell(1,nLags*nSources+1);
        tbltitles{1} = 't';
        
        n = 2;
        for s=1:nSources
            for l=1:nLags
                tbltitles{n} = ['n',num2str(sources(s)),'lag',num2str(abs(lags(l)))];
                
                n=n+1;
            end
        end
        
        X = array2table(X,'variablenames',tbltitles);
        X.Properties.VariableDescriptions{1} = 't';      
    end
        
    y = sinkRaster';
    y = y(:);
    
end