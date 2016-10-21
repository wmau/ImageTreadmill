function [edgesRemoved,newA] = ExciseOverlap(md,A)
%[edgesRemoved,newA] = ExciseOverlap(md,A)
%
%   Eliminates edges between overlapping neurons. 
%
%   INPUTS
%       md: session entry.
%
%       A: adjacency matrix.
%
%   OUPUTS
%       edgesRemoved: edge list of removed edges. 
%
%       newA: pruned adjacency matrix.
%

%% Set up.
    cd(md.Location); 
    load('FinalOutput.mat','NeuronPixels'); 
    
    el = adj2edgeL(A); el = el(:,1:2);
    nEdges = size(el,1);
    edgesRemoved = zeros(nEdges,2);
    
%% Remove edges from overlapping cells.
    newA = A;
    for e=1:length(el)
        src = el(e,1);
        snk = el(e,2); 
        
        %If there is even 1 pixel of overlap, remove. 
        if ~isempty(intersect(NeuronPixels{src},NeuronPixels{snk}))
            newA(src,snk) = 0;
            edgesRemoved(e,:) = el(e,:);
        end
    end
    
    %Clean up.
    edgesRemoved(edgesRemoved(:,1)==0,:) = [];
end