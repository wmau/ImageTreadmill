function [TCounts,TMap_gauss,TMap_unsmoothed] = ...
    MakePlacefield(spkpos,xEdges,yEdges,RunOccMap,varargin)
%[TCounts,TMap_gauss,TMap_unsmoothed] = ...
%    MakePlacefield(spkpos,xEdges,yEdges,RunOccMap,varargin)
%
%   Bite-sized place field function. Takes inputs after a bit of processing
%   to make place fields for one cell.
%
%   INPUTS
%       spkpos: 3xF (F = # frames) where the first row is FT from Tenaspis
%       output, the second row is x position, and the third row is y
%       position.
%   
%       xEdges & yEdges: output from MakeOccMap.
%
%       RunOccMap: output from MakeOccMap, occupancy map of where mouse
%       ran, in frame counts.
%
%       optional...
%           gauss_std: STD for gaussian filter (default = 2.5).
%
%           cmperbin: centimeters per spatial bin (default = 1). 
%
%   OUTPUTS
%       TCounts: Spike-place histogram (in frames counts).
%
%       TMap_gauss: Smoothed version of TMap_unsmoothed using Gaussian
%       filter.
%
%       TMap_unsmoothed: Spike-place histogram, normalized by RunOccMap.
%

%% Parse inputs.
    p = inputParser;
    p.addRequired('spkpos');
    p.addParameter('gauss_std',2.5,@(x) isscalar(x)); 
    p.addParameter('cmperbin',1,@(x) isscalar(x)); 
    
    p.parse(spkpos,varargin{:}); 
    
    gauss_std = p.Results.gauss_std; 
    cmperbin = p.Results.cmperbin;
    FT = logical(spkpos(1,:));
    x = spkpos(2,:);
    y = spkpos(3,:);
    
%% Make place field.
    gauss_std = gauss_std/cmperbin; 
    
    sm = fspecial('gaussian',[round(8*gauss_std,0),round(8*gauss_std)],gauss_std); 
    
    TCounts = hist2(y(FT),x(FT),yEdges,xEdges); 
    Tsum = sum(TCounts(:)); 
    
    %Normalize. 
    TMap_unsmoothed = TCounts./RunOccMap; 
    TMap_unsmoothed(isnan(TMap_unsmoothed)) = 0;
    
    %Smooth. 
    TMap_gauss = imfilter(TMap_unsmoothed,sm);
    TMap_gauss = TMap_gauss.*Tsum./sum(TMap_gauss(:)); 
    
    %Where the mouse didn't run, make nan. 
    TMap_unsmoothed(RunOccMap==0) = nan;
    TMap_gauss(RunOccMap==0) = nan;
end