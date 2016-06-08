function [LAGS,p,SHUFFLELAGS] = lapCC(raster1,raster2,B)
%LAGS = lapCC(md,n1,n2,shuffle)
%
%   Finds the "cross-correlation" between two neurons every treadmill run.
%   Specifically, it looks for onsets of calcium transients and calculates
%   differences between their timestamps (n1 - n2; so positive lags mean n2
%   occurred before n1). Then it fixes n1 in place and shuffles n2 in time,
%   computing the lag from this randomized data to create a null
%   distribution for statistical testing.
%
%   INPUTS
%       raster1 and raster2: Outputs from buildRaster for two neurons that
%       you want to compare. 
%
%   OUTPUTS
%       LAGS: Pairwise difference in seconds between each response in a lap
%       from n1 to n2. 

%% 
    [nLaps,nBins] = size(raster1);
    
%%    
    LAGS = [];
    SHUFFLELAGS = [];
    for l=1:nLaps
        if any(raster1(l,:)) && any(raster2(l,:))
            a = find(raster1(l,:))./20; 
            b = find(raster2(l,:))./20;

            lags = bsxfun(@minus,a,b');
            LAGS = [LAGS; lags(:)];

            for i=1:B
                Bb = find(circshift(raster2(l,:),[0,randi([0,nBins])]))./20;

                lags = bsxfun(@minus,a,Bb');
                SHUFFLELAGS = [SHUFFLELAGS; lags(:)];
            end
        end
    end
    
    if ~isempty(LAGS) && B > 0
        [~,p] = kstest2(LAGS,SHUFFLELAGS);
    else
        p = 1; 
    end
    
%     figure;
%     histogram(LAGS,[-200:5:200],'normalization','probability');
%     
%     if shuffle
%         hold on;
%         histogram(SHUFFLELAGS,[-200:5:200],'normalization','probability');
%         hold off;
%         [~,p] = kstest2(LAGS,SHUFFLELAGS);
%         title(num2str(p));
%     end

end