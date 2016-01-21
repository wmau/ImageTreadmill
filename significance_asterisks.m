function [SIGX,SIGY] = significance_asterisks(t,sigcurve,smoothedcurve,bins)
%[SIGX,SIGY] = significance_asterisks(t,sigcurve,smoothedcurve,bins)
%
%   Finds the X and Y locations of the smoothed tuning curve that is
%   statistically significant for time tuning. 
%
%   INPUTS
%       t: Time vector, same size as smoothed curve. 
%
%       sigcurve: From FindTimeCells. 
%
%       smoothedcurve: Smoothed tuning curve. 
%
%       bins: 1:nBins, same size as t. 
%
%   OUTPUTS
%       SIGX & SIGY: X and Y locations for where to place asterisk. 
%

%%
    %Get indices of significance.
    sigT = find(sigcurve);
    
    %Find the corresponding index in the smoothed bins.
    [~,inds] = ismember(sigT,bins);
    
    %Values to plot.
    SIGX = t(inds);
    SIGY = smoothedcurve(inds);
end