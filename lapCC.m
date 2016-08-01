function [LAGS,p,SHUFFLELAGS] = lapCC(trigRaster,targRaster,B)
%LAGS = lapCC(md,n1,n2,shuffle)
%
%   Finds the "cross-correlation" between two neurons every treadmill run.
%   Specifically, it looks for onsets of calcium transients and calculates
%   differences between their timestamps (n1 - n2; so positive lags mean n2
%   occurred before n1). Then it fixes n2 in place and shuffles n1 in time,
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
    [nLaps,nBins] = size(trigRaster);
    
%%    
    [~,LAGS] = stripRaster(trigRaster,targRaster);
    SHUFFLELAGS = [];
    
    if ~isempty(LAGS)
        for i=1:B
            temp1 = zeros(nLaps,nBins); 

            for l=1:nLaps
                temp1(l,:) = circshift(trigRaster(l,:),[0,randi([0,nBins])]);
            end

            [~,shuffledLags] = stripRaster(temp1,targRaster);
            SHUFFLELAGS = [SHUFFLELAGS; shuffledLags];
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