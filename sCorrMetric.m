function stability = sCorrMetric(normtilemat)
%stability = sCorrMetric(normtilemat)
%
%   Computes the Spearman correlation of each time cell in the first cell
%   array with itself on other cell arrays (other sessions). 
%
%   INPUT
%       normtilemat: output from multiPastalkovaPlot.
%
%   OUTPUT
%       stability: TCxS matrix of correlation coefficients for each time
%       cell across sessions. First column is always 1.  
%

%% Do correlations. 
    nSessions = length(normtilemat);
    nTCs = size(normtilemat{1},1);
    
    %Preallocate and assign first column all 1s. 
    stability = nan(nTCs,nSessions);
    stability(:,1) = 1; 
    for s=2:nSessions
        for n=1:nTCs
            %Spearman (non-parametric) correlation.
            stability(n,s) = corr(normtilemat{1}(n,:)',normtilemat{s}(n,:)','type','Pearson');
        end
    end

    %Set NaNs to 0. 
    stability(isnan(stability)) = 0;
end