function graphData = buildGraphData(md,varstring)
%
%
%

%% 
    cd(md.Location);
    load(fullfile(md.Location,varstring)); 
    
    graphData.A = A; 
    graphData.Ap = Ap; 
    graphData.CC = CC; 
    graphData.closest = closest;
    graphData.nulld = nulld;
    graphData.mdInfo = mdInfo;
    
end