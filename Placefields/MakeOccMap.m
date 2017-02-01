function [OccMap,RunOccMap,xEdges,yEdges,xBin,yBin] = ...
    MakeOccMap(x,y,lims,good,isrunning,cmperbin)
%[OccMap,RunOccMap,xEdges,yEdges] = ...
%    MakeOccMap(x,y,lims,good,isrunning,cmperbin)
%
%   Makes occupancy maps given X/Y limits and position. 
%
%   INPUTS
%       X & Y: mouse position, aligned.
%
%       lims: 2x2 matrix where [xmin xmax
%                               ymin ymax]
%
%       good: all frames minus ones specified by exclude frames. 
%
%       isrunning: all good frames where mouse velocity exceeds minspeed.
%
%       cmperbin: centimeteres per spatial bin.
%
%   OUPUTS
%       OccMap: occupancy map (in frames counts).
%
%       RunOccMap: occupancy map where mouse is running (in frame counts). 
%
%       xEdges & yEdges: edges used for histogram.
%

%% Extract limits.
    xmin = lims(1,1);
    xmax = lims(1,2); 
    ymin = lims(2,1);
    ymax = lims(2,2); 

%% Make edges for hist2.
    Xrange = xmax-xmin; 
    Yrange = ymax-ymin; 
    
    nXBins = ceil(Xrange/cmperbin); 
    nYBins = ceil(Yrange/cmperbin); 
    
    xEdges = (0:nXBins)*cmperbin+xmin;
    yEdges = (0:nYBins)*cmperbin+ymin; 
 
%% Run 2D histogram function.
    OccMap = histcounts2(x(good),y(good),xEdges,yEdges); 
    [RunOccMap,~,~,xBin,yBin] = histcounts2(x(isrunning),y(isrunning),xEdges,yEdges); 
end