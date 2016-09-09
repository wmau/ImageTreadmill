function saveGraphData(graphData,savename)
%saveGraphData(graphData,savename)
%
%   Saves graph variables separately in a mat file specified by savename. 
%
%   INPUTS
%       graphData: from MakeGraphv4 or pruneA.
%
%       savename: string, name of file to be saved.
%

%%
    md = findMDfromGraphData(graphData);
    cd(md.Location); 
    
    A = graphData.A; 
    Ap = graphData.Ap;
    CC = graphData.CC;
    nulld = graphData.nulld;
    closest = graphData.closest;
    mdInfo.Animal = graphData.Animal;
    mdInfo.Date = graphData.Date;
    mdInfo.Session = graphData.Session;
    
    if isfield(graphData,'prune_p');
        prune_p = graphData.prune_p; 
        save(savename,'A','Ap','nulld','closest','mdInfo','prune_p','-v7.3');
    else
        save(savename,'A','Ap','nulld','closest','mdInfo','-v7.3');      
    end
    
end