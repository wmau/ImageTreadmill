function [X,y] = SourceSinkGLMSetUp(md,sources,sink)
%
%
%

%% Set up.  
        cd(md.Location);
        load('Pos_align.mat','FT'); 
        load('TimeCells.mat','TimeCells','T','TodayTreadmillLog'); 

        %Make treadmill run indices even.
        inds = TodayTreadmillLog.inds; 
        inds = inds(find(TodayTreadmillLog.complete),:);
        inds(:,2) = inds(:,1) + 20*T-1; 

        sinkRaster = buildRaster(inds,FT,sink);
        [nLaps,nTimeBins] = size(sinkRaster); 

    if nargout > 1
        t = linspace(0,T,nTimeBins)';
        X(:,1) = repmat(t,nLaps,1);

        nSources = length(sources);
        for s=1:nSources
            sourceRaster = buildRaster(inds,FT,sources(s));
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
        
    y = sinkRaster(:);
    
end